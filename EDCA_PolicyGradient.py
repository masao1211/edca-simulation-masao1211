from numba import jit
import numpy as np
import numpy.random as rd
import copy
from functools import lru_cache
from joblib import Parallel, delayed
import itertools
import matplotlib.pyplot as plt
deepcopy = copy.deepcopy



# ////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////
# EDCA
# ////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////
def initialize_edca_magic_number(EDCAs):
    CWmin = np.array([3,7,15,15],dtype=int)
    CWmax = np.array([7,15,1023,1023],dtype=int)
    AIFSN = np.array([2,2,3,7],dtype=int)
    SIFS = 16
    SlotTime = 9
    SendTime = 248
    ACKTime = 24
    BusyTime = SIFS + SendTime + ACKTime
    return CWmin[:EDCAs], CWmax[:EDCAs], AIFSN[:EDCAs], SIFS, SlotTime, BusyTime
def do_edca(EDCAInput):
    PiFactor,ThetaStateActionElement,TimeSteps,EDCAs,APs,OtherAPs,ArrivalRate,InitialPacketInQueue = EDCAInput
    EDCAMagicNumber = initialize_edca_magic_number(EDCAs)
    CWmin, CWmax, AIFSN, SIFS, SlotTime, BusyTime = EDCAMagicNumber
    QueueState,CountState,OtherCountState,ToSendPacket,Avo,Nvo,time,StateActionMemo,NablaLogPis = initialize_edca(EDCAInput,EDCAMagicNumber)
    # print(time,Avo,Nvo,QueueState[1],CountState,ToSendPacket)
    # print(time,Avo,Nvo,QueueState[1],CountState,ToSendPacket,OtherCountState)
    for t in range(10000):
        QueueState,CountState,OtherCountState,ToSendPacket,Avo,Nvo,time,StateActionMemo,NablaLogPis = arrrive_and_wait_and_send(EDCAInput,EDCAMagicNumber,QueueState,CountState,OtherCountState,ToSendPacket,Avo,Nvo,time,StateActionMemo,NablaLogPis)
        # print(time,Avo,Nvo,QueueState[1],CountState,ToSendPacket)
        # print(time,Avo,Nvo,QueueState[1],CountState,ToSendPacket,OtherCountState)
        if np.all(Nvo >= TimeSteps): return time-SIFS,StateActionMemo,np.sum(NablaLogPis,axis=1)
    print(0)
    return 0
def initialize_edca(EDCAInput,EDCAMagicNumber):
    PiFactor,ThetaStateActionElement,TimeSteps,EDCAs,APs,OtherAPs,ArrivalRate,InitialPacketInQueue = EDCAInput
    CWmin, CWmax, AIFSN, SIFS, SlotTime, BusyTime = EDCAMagicNumber
    ToSendPacket = np.ones((APs,EDCAs),dtype=int)
    PacketQueue = np.array([deepcopy(initialize_packet_queue(TimeSteps,tuple(InitialPacketInQueue))) for ap in range(APs)],dtype=int)
    PacketsInQueue = np.array([InitialPacketInQueue for ap in range(APs)],dtype=int)
    CW = np.array([deepcopy(CWmin) for ap in range(APs)],dtype=int)
    Avo = np.zeros(APs,dtype=int)
    StateActionMemoElement = (APs,TimeSteps,APs*5+OtherAPs*2+1)
    StateActionMemo = np.empty(StateActionMemoElement,dtype=int)
    NablaLogPis = np.zeros((APs,TimeSteps) + ThetaStateActionElement,dtype=float)

    BackoffCount = np.zeros((APs,EDCAs),dtype=int)
    OtherBackoffCount = np.zeros((OtherAPs,EDCAs),dtype=int)

    ArrivePacket = np.zeros((APs,EDCAs),dtype=int)
    ArrivePacket[0][0] = 1
    QueueState = PacketQueue,PacketsInQueue
    QueueState,Avo,StateActionMemo,NablaLogPis = ArrivePacket_to_PacketQueue(EDCAInput,ArrivePacket,QueueState,BackoffCount,OtherBackoffCount,Avo,StateActionMemo,NablaLogPis)

    AdvanceTime = SIFS
    QueueState,Avo,StateActionMemo,NablaLogPis = arrive_packet(EDCAInput,QueueState,BackoffCount,OtherBackoffCount,AdvanceTime,Avo,StateActionMemo,NablaLogPis)

    AIFSCount = np.array([[aifsn if piq>0 else -1 for piq,aifsn in zip(PacketsInQueue[ap],AIFSN)] for ap in range(APs)],dtype=int)
    BackoffCount = np.array([[rd.randint(cw+1) if piq>0 else 0 for piq,cw in zip(PacketsInQueue[ap],CW[ap])] for ap in range(APs)],dtype=int)

    OtherCW = np.array([deepcopy(CWmin) for i in range(OtherAPs)],dtype=int)
    OtherAIFSCount = np.array([deepcopy(AIFSN) for i in range(OtherAPs)],dtype=int)
    OtherBackoffCount = np.array([[rd.randint(cw+1) for cw in cw1] for cw1 in OtherCW],dtype=int)

    Nvo = np.zeros(APs,dtype=int)
    time = SIFS
    CountState = BackoffCount,AIFSCount,CW
    OtherCountState = OtherBackoffCount,OtherAIFSCount,OtherCW

    return QueueState,CountState,OtherCountState,ToSendPacket,Avo,Nvo,time,StateActionMemo,NablaLogPis
# @lru_cache(maxsize=1)
def initialize_packet_queue(TimeSteps,InitialPacketInQueue):
    PacketQueue = np.array([[(j <= iq)*5-1 for j in range(TimeSteps*10)] for iq in InitialPacketInQueue],dtype=int)
    return PacketQueue

# ////////////////////////////////////////////////////////////////////////////////////////////////////
# arrive packets
# ////////////////////////////////////////////////////////////////////////////////////////////////////
def arrive_packet(EDCAInput,QueueState,BackoffCount,OtherBackoffCount,AdvanceTime,Avo,StateActionMemo,NablaLogPis):
    PiFactor,ThetaStateActionElement,TimeSteps,EDCAs,APs,OtherAPs,ArrivalRate,InitialPacketInQueue = EDCAInput
    ArrivePacket = np.array([[rd.poisson(ar*AdvanceTime) for ar in ArrivalRate] for ap in range(APs)],dtype=int)
    QueueState,Avo,StateActionMemo,NablaLogPis = ArrivePacket_to_PacketQueue(EDCAInput,ArrivePacket,QueueState,BackoffCount,OtherBackoffCount,Avo,StateActionMemo,NablaLogPis)
    return QueueState,Avo,StateActionMemo,NablaLogPis
def ArrivePacket_to_PacketQueue(EDCAInput,ArrivePacket,QueueState,BackoffCount,OtherBackoffCount,Avo,StateActionMemo,NablaLogPis):
    PiFactor,ThetaStateActionElement,TimeSteps,EDCAs,APs,OtherAPs,ArrivalRate,InitialPacketInQueue = EDCAInput
    PacketQueue,PacketsInQueue = QueueState
    Array = [[[] for i in range(EDCAs)] for ap in range(APs)]
    for i in range(np.sum(ArrivePacket)):
        APChoice,PacketChoice = choose_queue(EDCAs,APs,ArrivePacket)
        ArrivePacket[APChoice][PacketChoice] -= 1
        if PacketChoice == 0 and Avo[APChoice] < TimeSteps:
            PacketsInQueueList = [piq for piq1 in PacketsInQueue for piq in piq1]
            BackoffCountList = [bc for bc1 in BackoffCount for bc in bc1]
            OtherBackoffCountList = [obc for obc1 in OtherBackoffCount for obc in obc1]
            StateList = PacketsInQueueList + BackoffCountList + OtherBackoffCountList + Avo.tolist()
            state = tuple(StateList)
            Phi = calculate_phi_no_indicator(state,PiFactor,ThetaStateActionElement)
            Pi,NablaLogPiSecondMember = calculate_pi_and_second_member(state,PiFactor,ThetaStateActionElement,EDCAs,Phi,APChoice)
            QueueChoice = rd.multinomial(1, Pi)
            QueueChoice = np.sum([i*c for i,c in zip(range(EDCAs),QueueChoice)])

            state_action_memo = np.array(StateList + [QueueChoice],dtype=int)
            StateActionMemo[APChoice][Avo[APChoice]] = state_action_memo
            NablaLogPis[APChoice][Avo[APChoice]] = calculate_nabla_log_pi(state,QueueChoice,PiFactor,ThetaStateActionElement,EDCAs,Phi,NablaLogPiSecondMember)

            Array[APChoice][QueueChoice].append(0)
            PacketsInQueue[APChoice][QueueChoice] += 1
            Avo[APChoice] += 1
        else:
            Array[APChoice][PacketChoice].append(PacketChoice)
            PacketsInQueue[APChoice][PacketChoice] += 1
    AddPacketsInQueue = np.array([[len(a) for a in Array[ap]] for ap in range(APs)],dtype=int)
    Nmin = np.array([[np.where(pq==-1)[0][0] for pq in PacketQueue[ap]] for ap in range(APs)],dtype=int)
    Nmax = np.array([[n+apiq if n+apiq<TimeSteps*10-1 else TimeSteps*10-1 for n,apiq in zip(Nmin[ap],AddPacketsInQueue[ap])] for ap in range(APs)],dtype=int)
    PacketQueue = np.array([get_packet_queue(PacketQueue[ap],Nmin[ap],Nmax[ap],Array[ap],EDCAs) for ap in range(APs)],dtype=int)
    QueueState = PacketQueue,PacketsInQueue
    return QueueState,Avo,StateActionMemo,NablaLogPis

def calculate_phi_no_indicator(state,PiFactor,ThetaStateActionElement):
    Theta,ThetaStateElement,ThetaStateIterator,PhiElementParams = PiFactor
    # Phi = Phi(state) state:固定,action:指示関数が0にならないようなaction
    Phi = np.zeros(ThetaStateElement,dtype=float)
    state_array = np.array(state,dtype=float)
    PhiElementAdd,PhiElementProd = PhiElementParams
    for tsi in ThetaStateIterator:
        Phi[tsi] = np.prod((PhiElementProd*(state_array+PhiElementAdd))**np.array(tsi,dtype=float))
    return Phi

def calculate_pi_and_second_member(state,PiFactor,ThetaStateActionElement,EDCAs,Phi,APChoice):
    Theta,ThetaStateElement,ThetaStateIterator,PhiElementParams = PiFactor
    # ThetaPhi = [ThetaPhi(state,VO),ThetaPhi(state,VI)] state:固定,actionの数の配列
    ThetaPhi = np.array([np.sum([Theta[APChoice][tsi][action]*Phi[tsi] for tsi in ThetaStateIterator]) for action in range(EDCAs)],dtype=float)
    ExpThetaPhi = np.exp(ThetaPhi-np.max(ThetaPhi))
    SumExpThetaPhi = np.sum(ExpThetaPhi)
    Pi = np.array([ExpThetaPhi[action]/SumExpThetaPhi for action in range(EDCAs)],dtype=float)
    # //////////////////////////////////////////////////////////////////////////
    SumPhiExpThetaPhi = np.zeros(ThetaStateActionElement,dtype=float)
    for tsi in ThetaStateIterator:
        SumPhiExpThetaPhi[tsi] = np.array([np.sum([Phi[tsi]*ExpThetaPhi[action] if LastSub==action else 0 for action in range(EDCAs)]) for LastSub in range(EDCAs)],dtype=float)
    NablaLogPiSecondMember = SumPhiExpThetaPhi / SumExpThetaPhi
    return Pi,NablaLogPiSecondMember


def calculate_nabla_log_pi(state,QueueChoice,PiFactor,ThetaStateActionElement,EDCAs,Phi,NablaLogPiSecondMember):
    Theta,ThetaStateElement,ThetaStateIterator,PhiElementParams = PiFactor
    NablaLogPi = np.zeros(ThetaStateActionElement,dtype=float)
    for tsi in ThetaStateIterator:
        NablaLogPi[tsi] = np.array([Phi[tsi] - NablaLogPiSecondMember[tsi][LastSub] if QueueChoice==LastSub else -NablaLogPiSecondMember[tsi][LastSub] for LastSub in range(EDCAs)],dtype=float)
    return NablaLogPi

# print(1-np.ones(10))

def choose_queue(EDCAs,APs,ArrivePacket):
    r = rd.randint(np.sum(ArrivePacket))
    for i in range(APs):
        for j in range(EDCAs):
            r -= ArrivePacket[i][j]
            if r<0: return i,j
def get_packet_queue(PacketQueue,Nmin,Nmax,Array,EDCAs):
    for i in range(EDCAs):
        PacketQueue[i][Nmin[i]:Nmax[i]] = Array[i][:Nmax[i]-Nmin[i]]
    return PacketQueue
def append_list(array, EDCAs):
    a = [j * np.ones(i,dtype=int) for i,j in zip(array,range(EDCAs))]
    b = [ i for arr in a for i in arr]
    return b

# ////////////////////////////////////////////////////////////////////////////////////////////////////
# reduce counts and send
# ////////////////////////////////////////////////////////////////////////////////////////////////////
def arrrive_and_wait_and_send(EDCAInput,EDCAMagicNumber,QueueState,CountState,OtherCountState,ToSendPacket,Avo,Nvo,time,StateActionMemo,NablaLogPis):
    PiFactor,ThetaStateActionElement,TimeSteps,EDCAs,APs,OtherAPs,ArrivalRate,InitialPacketInQueue = EDCAInput
    CountState,OtherCountState,AdvanceTime,Avo,StateActionMemo,NablaLogPis = reduce_count(EDCAInput,EDCAMagicNumber,QueueState,CountState,OtherCountState,Avo,StateActionMemo,NablaLogPis)
    BackoffCount,AIFSCount,CW = CountState
    OtherBackoffCount,OtherAIFSCount,OtherCW = OtherCountState
    if np.any(AIFSCount**2+BackoffCount**2==0) or np.any(OtherAIFSCount**2+OtherBackoffCount**2==0):
        QueueState,CountState,OtherCountState,ToSendPacket,Avo,Nvo,AdvanceTime,StateActionMemo,NablaLogPis = do_send(EDCAInput,EDCAMagicNumber,QueueState,CountState,OtherCountState,ToSendPacket,Avo,Nvo,AdvanceTime,StateActionMemo,NablaLogPis)
    time += AdvanceTime
    return QueueState,CountState,OtherCountState,ToSendPacket,Avo,Nvo,time,StateActionMemo,NablaLogPis
def reduce_count(EDCAInput,EDCAMagicNumber,QueueState,CountState,OtherCountState,Avo,StateActionMemo,NablaLogPis):
    PiFactor,ThetaStateActionElement,TimeSteps,EDCAs,APs,OtherAPs,ArrivalRate,InitialPacketInQueue = EDCAInput
    CWmin, CWmax, AIFSN, SIFS, SlotTime, BusyTime = EDCAMagicNumber
    BackoffCount,AIFSCount,CW = CountState
    OtherBackoffCount,OtherAIFSCount,OtherCW = OtherCountState
    if np.any(Avo<TimeSteps):
        AdvanceTime = SlotTime
        QueueState,Avo,StateActionMemo,NablaLogPis = arrive_packet(EDCAInput,QueueState,BackoffCount,OtherBackoffCount,AdvanceTime,Avo,StateActionMemo,NablaLogPis)
        BackoffCount = np.array([[bc-1 if ac==0 else bc for ac,bc in zip(AIFSCount[ap],BackoffCount[ap])] for ap in range(APs)],dtype=int)
        AIFSCount = np.array([[ac-1 if ac>0 else ac for ac in AIFSCount[ap]] for ap in range(APs)],dtype=int)
        OtherBackoffCount = np.array([[bc-1 if ac==0 else bc for ac,bc in zip(ac1,bc1)] for ac1,bc1 in zip(OtherAIFSCount,OtherBackoffCount)],dtype=int)
        OtherAIFSCount = np.array([[ac-1 if ac>0 else ac for ac in ac1] for ac1 in OtherAIFSCount],dtype=int)
        if OtherAPs==0 and np.all(AIFSCount==-np.ones((APs,EDCAs),dtype=int)):
            PacketQueue,PacketsInQueue = QueueState
            AIFSCount = np.array([[aifsn if piq>0 else -1  for piq,aifsn in zip(PacketsInQueue[ap],AIFSN)] for ap in range(APs)],dtype=int)
            BackoffCount = np.array([[0 if piq==0 else (rd.randint(cw+1) if bc==0 else bc) for piq,cw,bc in zip(PacketsInQueue[ap],CW[ap],BackoffCount[ap])] for ap in range(APs)],dtype=int)
    else:
        AIFSPlusBackoffCountMin = np.min([[10000 if ac==-1 else ac+bc for ac,bc in zip(AIFSCount[ap],BackoffCount[ap])] for ap in range(APs)])
        OtherAIFSPlusBackoffCount = np.min(OtherBackoffCount+OtherAIFSCount) if OtherAPs>0 else 10000
        ReduceCount = AIFSPlusBackoffCountMin if AIFSPlusBackoffCountMin<OtherAIFSPlusBackoffCount else OtherAIFSPlusBackoffCount
        AdvanceTime = SlotTime * ReduceCount
        QueueState,Avo,StateActionMemo,NablaLogPis = arrive_packet(EDCAInput,QueueState,BackoffCount,OtherBackoffCount,AdvanceTime,Avo,StateActionMemo,NablaLogPis)
        BackoffCount = np.array([[0 if ac==-1 else bc-(ReduceCount-ac)*(ReduceCount-ac>0) for ac,bc in zip(AIFSCount[ap],BackoffCount[ap])] for ap in range(APs)],dtype=int)
        AIFSCount = np.array([[-1 if ac==-1 else (ac-ReduceCount)*(ac-ReduceCount>0) for ac in AIFSCount[ap]] for ap in range(APs)],dtype=int)
        OtherBackoffCount = np.array([[bc-(ReduceCount-ac)*(ReduceCount-ac>0) for ac,bc in zip(ac1,bc1)] for ac1,bc1 in zip(OtherAIFSCount,OtherBackoffCount)],dtype=int)
        OtherAIFSCount = np.array([[(ac-ReduceCount)*(ac-ReduceCount>0) for ac,bc in zip(ac1,bc1)] for ac1,bc1 in zip(OtherAIFSCount,OtherBackoffCount)],dtype=int)

    CountState = BackoffCount,AIFSCount,CW
    OtherCountState = OtherBackoffCount,OtherAIFSCount,OtherCW
    return CountState,OtherCountState,AdvanceTime,Avo,StateActionMemo,NablaLogPis
def do_send(EDCAInput,EDCAMagicNumber,QueueState,CountState,OtherCountState,ToSendPacket,Avo,Nvo,AdvanceTime,StateActionMemo,NablaLogPis):
    PiFactor,ThetaStateActionElement,TimeSteps,EDCAs,APs,OtherAPs,ArrivalRate,InitialPacketInQueue = EDCAInput
    PacketQueue,PacketsInQueue = QueueState
    BackoffCount,AIFSCount,CW = CountState
    OtherBackoffCount,OtherAIFSCount,OtherCW = OtherCountState
    CWmin,CWmax,AIFSN,SIFS,SlotTime,BusyTime = EDCAMagicNumber
    SendQueue = BackoffCount**2+AIFSCount**2==np.zeros((APs,EDCAs),dtype=int)
    SendSTA = np.array([np.any(sq) for sq in SendQueue],dtype=int)
    if OtherAPs>0:
        OtherSendQueue = OtherBackoffCount**2+OtherAIFSCount**2==np.zeros((OtherAPs,EDCAs),dtype=int)
        OtherSendSTA = np.array([np.any(osq) for osq in OtherSendQueue],dtype=int)
    else:
        OtherSendSTA = 0
    SumSendSTA = np.sum(SendSTA)
    SumSendOtherSTA = np.sum(OtherSendSTA)
    collision = SumSendSTA + SumSendOtherSTA > 1

    if not collision:
        if SumSendSTA:
            SendSTANumber = np.where(SendSTA)[0][0]
            CanSendQueueNumber = np.where(SendQueue[SendSTANumber])[0][0]
            Nvo[SendSTANumber] = Nvo[SendSTANumber] + 1 if PacketQueue[SendSTANumber][CanSendQueueNumber][ToSendPacket[SendSTANumber][CanSendQueueNumber]] == 0 else Nvo[SendSTANumber]
            PacketsInQueue[SendSTANumber][CanSendQueueNumber] -= 1
            ToSendPacket[SendSTANumber][CanSendQueueNumber] += 1
            SendQueue[SendSTANumber][CanSendQueueNumber] = 0
            CW[SendSTANumber][CanSendQueueNumber] = CWmin[CanSendQueueNumber]
        else:
            SendSTANumber = np.where(OtherSendSTA)[0][0]
            CanSendQueueNumber = np.where(OtherSendQueue[SendSTANumber])[0][0]
            OtherSendQueue[SendSTANumber][CanSendQueueNumber] = 0
            OtherCW[SendSTANumber][CanSendQueueNumber] = CWmin[CanSendQueueNumber]

    AdvanceTime += BusyTime + SIFS
    if np.any(Avo<TimeSteps):
        QueueState = PacketQueue,PacketsInQueue
        QueueState,Avo,StateActionMemo,NablaLogPis = arrive_packet(EDCAInput,QueueState,BackoffCount,OtherBackoffCount,AdvanceTime,Avo,StateActionMemo,NablaLogPis)
        PacketQueue,PacketsInQueue = QueueState

    CW = np.array([[(2*cw+1 if 2*cw+1<cwmax else cwmax) if sq else cw for cw,cwmax,sq in zip(CW[ap],CWmax,SendQueue[ap])] for ap in range(APs)],dtype=int)
    AIFSCount = np.array([[aifsn if piq>0 else -1  for piq,aifsn in zip(PacketsInQueue[ap],AIFSN)] for ap in range(APs)],dtype=int)
    BackoffCount = np.array([[0 if piq==0 else (rd.randint(cw+1) if bc==0 else bc) for piq,cw,bc in zip(PacketsInQueue[ap],CW[ap],BackoffCount[ap])] for ap in range(APs)],dtype=int)

    if OtherAPs>0:
        OtherCW = np.array([[(2*cw+1 if 2*cw+1<cwmax else cwmax) if sq else cw for cw,cwmax,sq in zip(cw1,CWmax,sq1)] for cw1,sq1 in zip(OtherCW,OtherSendQueue)],dtype=int)
        OtherAIFSCount = np.array([deepcopy(AIFSN) for i in range(OtherAPs)],dtype=int)
        OtherBackoffCount = np.array([[rd.randint(cw+1) if bc==0 else bc for cw,bc in zip(cw1,bc1)] for cw1,bc1 in zip(OtherCW,OtherBackoffCount)],dtype=int)

    QueueState = PacketQueue,PacketsInQueue
    CountState = BackoffCount,AIFSCount,CW
    OtherCountState = OtherBackoffCount,OtherAIFSCount,OtherCW
    return QueueState,CountState,OtherCountState,ToSendPacket,Avo,Nvo,AdvanceTime,StateActionMemo,NablaLogPis

# ////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////
# static policies
# ////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////
def delay_of_static_policy(Repeat,mode,EDCAs,APs,OtherAPs,TimeSteps,ArrivalRate,InitialPacketInQueue):
    # rd.seed(seed=1)
    EDCAInput = initialize_edca_input(mode,EDCAs,APs,OtherAPs,TimeSteps,tuple(ArrivalRate),tuple(InitialPacketInQueue))
    delay,StateActionMemo,SumNablaLogPi = get_delay_etc(EDCAInput,Repeat)
    return np.sum(delay) / Repeat
    # return np.sum(delay) / Repeat, StateActionMemo
    # return np.sum(delay) / Repeat, StateActionMemo, SumNablaLogPi
def initialize_edca_input(mode,EDCAs,APs,OtherAPs,TimeSteps,ArrivalRate,InitialPacketInQueue):
    Degrees = 2
    ThetaStateActionElement,ThetaStateIterator,ThetaStateElement = calculate_theta_element_iterator(Degrees,EDCAs,APs,OtherAPs)
    Theta = np.zeros((APs,) + ThetaStateActionElement,dtype=float)
    PhiElementParams = (np.zeros(APs*5+OtherAPs*2,dtype=int),np.ones(APs*5+OtherAPs*2,dtype=int))
    if mode=='VO':
        # 必ずVO
        for ap in range(APs):
            Theta[ap][ThetaStateIterator[0]] = np.array([1000,-1000])
    elif mode=='VI':
        # 必ずVI
        for ap in range(APs):
            Theta[ap][ThetaStateIterator[0]] = np.array([-1000,1000])
    elif mode=='1/n':
        # 確率0.5
        for ap in range(APs):
            Theta[ap][ThetaStateIterator[0]] = np.array([0,0])
    elif mode=='less':
        # 少ない方
        for ap in range(APs):
            Theta[ap][ThetaStateIterator[0]] = np.array([100,-100])
            Theta[ap][ThetaStateIterator[ap*2+1]] = np.array([-1000000,1000000])
            Theta[ap][ThetaStateIterator[ap*2+2]] = np.array([1010000,-1010000])
    ArrivalRate = ArrivalRate[:EDCAs]
    InitialPacketInQueue = InitialPacketInQueue[:EDCAs]
    PiFactor = Theta,ThetaStateElement,ThetaStateIterator,PhiElementParams
    return PiFactor,ThetaStateActionElement,TimeSteps,EDCAs,APs,OtherAPs,ArrivalRate,InitialPacketInQueue
def calculate_theta_element_iterator(Degrees,EDCAs,APs,OtherAPs):
    ThetaChoice = [itertools.combinations_with_replacement(range(APs*5+OtherAPs*2), d) for d in range(Degrees)]
    ThetaChoiceOneHot = [[np.identity(APs*5+OtherAPs*2,dtype=int)[list(ti)] for ti in tc] for tc in ThetaChoice]
    ThetaStateIterators = [[tuple(np.sum(tioh,axis=0)) for tioh in tcoh] for tcoh in ThetaChoiceOneHot]
    ThetaStateIterator = list(itertools.chain.from_iterable(ThetaStateIterators))
    ThetaStateActionElement = (Degrees,)*(APs*5+OtherAPs*2) + (EDCAs,)
    ThetaStateElement = (Degrees,)*(APs*5+OtherAPs*2)
    return ThetaStateActionElement,ThetaStateIterator,ThetaStateElement
def get_delay_etc(EDCAInput,Repeat):
    if Repeat>=100:
        # delay,StateActionMemo,SumNablaLogPi = zip(*Parallel(n_jobs=-1)( [delayed(do_edca)(EDCAInput) for i in range(Repeat)] ))
        delay,StateActionMemo,SumNablaLogPi = zip(*Parallel(n_jobs=-1)( [delayed(do_edca_sub)(EDCAInput) for i in range(Repeat//100)] ))
        delay = np.array([delay3 for delay2 in delay for delay3 in delay2],dtype=int)
        StateActionMemo = np.array([StateActionMemo3 for StateActionMemo2 in StateActionMemo for StateActionMemo3 in StateActionMemo2],dtype=int)
        SumNablaLogPi = np.array([SumNablaLogPi3 for SumNablaLogPi2 in SumNablaLogPi for SumNablaLogPi3 in SumNablaLogPi2],dtype=float)
    else:
        delay,StateActionMemo,SumNablaLogPi = zip(*[do_edca(EDCAInput) for i in range(Repeat)])
    return delay,StateActionMemo,SumNablaLogPi
def do_edca_sub(EDCAInput):
    delay,StateActionMemo,SumNablaLogPi = zip(*[do_edca(EDCAInput) for i in range(100)])
    # delay,StateActionMemo,SumNablaLogPi = zip(*Parallel(n_jobs=-1)( [delayed(do_edca)(EDCAInput) for i in range(100)] ))
    return delay,StateActionMemo,SumNablaLogPi


# ThetaStateActionElement,ThetaStateIterator,ThetaStateElement = calculate_theta_element_iterator(Degrees=2,EDCAs=2,APs=2,OtherAPs=0)
# print(ThetaStateIterator)
# print(ThetaStateElement)


# delay,StateActionMemo,SumNablaLogPi = zip(*[do_edca(EDCAInput) for i in range(100)])
# delay,StateActionMemo,SumNablaLogPi = do_many_edca(EDCAInput)
# delay,StateActionMemo,SumNablaLogPi = zip(*Parallel(n_jobs=-1)( [delayed(do_many_edca)(EDCAInput) for i in range(2)] ))

# Repeat = 1000
# print(np.array(delay))

# 方策勾配

# ////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////
# 方策勾配法
# ////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////
def policy_gradient(Updates,Episodes,TimeSteps,Degrees,EDCAs,APs,OtherAPs,LearningRate,ArrivalRate,InitialPacketInQueue,PhiElementParams):
    ThetaStateActionElement,ThetaStateIterator,ThetaStateElement = calculate_theta_element_iterator(Degrees,EDCAs,APs,OtherAPs)
    # Thetas = np.zeros((Updates+1,APs) + ThetaStateActionElement,dtype=float)
    ThetaNow = np.zeros((APs,) + ThetaStateActionElement,dtype=float)
    ThetaNext = np.zeros((APs,) + ThetaStateActionElement,dtype=float)
    JThetas = np.empty(Updates+1,dtype=int)
    ArrivalRate = ArrivalRate[:EDCAs]
    InitialPacketInQueue = InitialPacketInQueue[:EDCAs]
    for k in range(Updates+1):
        ThetaNow = ThetaNext
        # PiFactor = Thetas[k],ThetaStateElement,ThetaStateIterator,PhiElementParams
        PiFactor = ThetaNow,ThetaStateElement,ThetaStateIterator,PhiElementParams
        EDCAInput = PiFactor,ThetaStateActionElement,TimeSteps,EDCAs,APs,OtherAPs,ArrivalRate,InitialPacketInQueue
        delays,StateActionMemos,SumNablaLogPis = get_delay_etc(EDCAInput,Episodes)
        JThetas[k] = np.average(delays)
        if k == Updates: break
        # BHatStar = Parallel(n_jobs=-1)( [delayed(calculate_b_hat_star)(delays,SumNablaLogPis,Episodes,ap) for ap in range(APs)] )
        BHatStar = np.array([calculate_b_hat_star(delays,SumNablaLogPis,Episodes,ap) for ap in range(APs)],dtype=float)
        BaselineDelays = np.array([delays - BHatStar[ap] for ap in range(APs)],dtype=float)
        # NablaJTheta = Parallel(n_jobs=-1)( [delayed(calculate_b_hat_star)(delays,SumNablaLogPis,Episodes,ap) for ap in range(APs)] )
        NablaJTheta = np.array([np.average([BaselineDelays[ap][m]*SumNablaLogPis[m][ap] for m in range(Episodes)],axis=0) for ap in range(APs)],dtype=float)
        # Thetas[k+1] = Thetas[k] - LearningRate * NablaJTheta
        ThetaNext = ThetaNow - LearningRate * NablaJTheta
    # return JThetas, Thetas
    return JThetas, ThetaNow
    # return JThetas,
def calculate_b_hat_star(delays,SumNablaLogPis,Episodes,ap):
    numerator = np.sum([delays[m]*SumNablaLogPis[m][ap]*SumNablaLogPis[m][ap] for m in range(Episodes)])
    denominator = np.sum([SumNablaLogPis[m][ap]*SumNablaLogPis[m][ap] for m in range(Episodes)])
    return numerator/denominator


# ////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////
# 方策の計算
# ////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////
def calculate_pi(Theta,state,PiFactor,ThetaStateActionElement,ThetaStateIterator,EDCAs,APChoice):
    Phi = calculate_phi_no_indicator(state,PiFactor,ThetaStateActionElement)
    ThetaPhi = np.array([np.sum([Theta[APChoice][tsi][action]*Phi[tsi] for tsi in ThetaStateIterator]) for action in range(EDCAs)],dtype=float)
    ExpThetaPhi = np.exp(ThetaPhi-np.max(ThetaPhi))
    SumExpThetaPhi = np.sum(ExpThetaPhi)
    Pi = np.array([ExpThetaPhi[action]/SumExpThetaPhi for action in range(EDCAs)],dtype=float)
    return Pi


# ////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////
# 方策の図示
# ////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////
def calculate_and_save_result(Updates,Episodes,TimeSteps,Degrees,EDCAs,APs,OtherAPs,LearningRate,ArrivalRate,InitialPacketInQueue,PhiElementParams):
    JThetas, Thetas = policy_gradient(Updates,Episodes,TimeSteps,Degrees,EDCAs,APs,OtherAPs,LearningRate,ArrivalRate,InitialPacketInQueue,PhiElementParams)
    np.save('JThetas',JThetas)
    np.save('Thetas',Thetas)
def load_and_show_learning_curve(Updates,APs,OtherAPs):
    JThetas = np.load('JThetas.npy')
    JThetas = JThetas[0:Updates+1]
    x = np.arange(Updates + 1)
    VODelay = [[0,0,0],[3895,11495,20807],[12651,0,0]]
    HeuristicsDelay = [[0,0,0],[3804,10424,18298],[11507,0,0]]
    PiVO = VODelay[APs][OtherAPs] * np.ones(Updates + 1)
    PiHeuristics = HeuristicsDelay[APs][OtherAPs] * np.ones(Updates + 1)
    # plt.plot(x, JThetas)
    plt.plot(x, JThetas, label="Proposed")
    plt.plot(x, PiVO, label="Conventional")
    plt.plot(x, PiHeuristics, label="Heuristics")
    plt.xlim([0,Updates])
    # plt.ylim([10000,20000])
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.xlabel('Number of updating')
    # plt.xlabel('Number of updating θ')
    plt.ylabel("Transmission delay time (μs)")
    plt.legend(loc='upper right', borderaxespad=2)
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    # plt.subplots_adjust(right=0.66)
    plt.savefig('learning_curve.pdf')
def load_and_show_policy(Updates,Degrees,EDCAs,APs,OtherAPs,PhiElementParams):
    ThetaStateActionElement,ThetaStateIterator,ThetaStateElement = calculate_theta_element_iterator(Degrees,EDCAs,APs,OtherAPs)

    Theta = np.load('Thetas.npy')
    # Thetas = np.load('Thetas.npy')
    # Theta = Thetas[Updates]

    # Degrees = 2
    # ThetaStateActionElement,ThetaStateIterator,ThetaStateElement = calculate_theta_element_iterator(Degrees,EDCAs,APs,OtherAPs)
    # Theta = np.zeros((APs,) + ThetaStateActionElement,dtype=float)
    # PhiElementParams = (np.zeros(APs*5+OtherAPs*2,dtype=int),np.ones(APs*5+OtherAPs*2,dtype=int))
    # for ap in range(APs):
    #     Theta[ap][ThetaStateIterator[0]] = np.array([100,-100])
    #     Theta[ap][ThetaStateIterator[ap*2+1]] = np.array([-1000000,1000000])
    #     Theta[ap][ThetaStateIterator[ap*2+2]] = np.array([1010000,-1010000])

    PiFactor = Theta,ThetaStateElement,ThetaStateIterator,PhiElementParams

    x = np.arange(10)
    for Qvo in [0,4,8]:
        y = np.empty(10)
        for Qvi in range(10):
            # state = [0,0,0,0,0]*APs + [0,0]*OtherAPs
            state = [3,3,3,3,3,3,3,3,3,3]
            state[8] = Qvo
            state[9] = Qvi
            state = tuple(state)
            APChoice = 0
            Pi = calculate_pi(Theta,state,PiFactor,ThetaStateActionElement,ThetaStateIterator,EDCAs,APChoice)
            y[Qvi] = Pi[0]

        # plt.plot(x, y, label="$S_2 =$" + "{}".format(Qvo))
        # plt.plot(x, y, label="$C_\mathrm{VO}^{(1)} =$" + "{}".format(Qvo))
        plt.plot(x, y, label="$S_1 =$" + "{}".format(Qvo))

    plt.xlim([0,9])
    plt.ylim([0,1.01])

    # plt.xlabel("Number of packets in VI queue of AP 1")
    # plt.xlabel("Number of backoff count in VI queue of AP 1")
    plt.xlabel("Number of packet having arrived at AP 2")

    plt.ylabel("Probability of mapping to VO queue of AP 1")
    # plt.ylabel("Probability of mapping to VO queue of AP 2")

    plt.legend(loc='upper left', borderaxespad=1)

    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    # plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0,title="$Q_\mathrm{VO}$ shows\nnumber\nof packets\nin VO queue")
    # plt.subplots_adjust(right=0.77)
    plt.savefig('figure-VO.pdf')
def delay_of_gotten_policy(Repeat,Updates,Degrees,PhiElementParams,EDCAs,OtherAPs,TimeSteps,ArrivalRate,InitialPacketInQueue):
    rd.seed(seed=1)
    delay = 0
    Thetas = np.load('Thetas.npy')
    Theta = Thetas[Updates]
    ArrivalRate = ArrivalRate[:EDCAs]
    InitialPacketInQueue = InitialPacketInQueue[:EDCAs]
    ThetaStateActionElement,ThetaStateIterator,ThetaStateElement = calculate_theta_element_iterator(Degrees,EDCAs,OtherAPs)
    PiFactor = Theta,ThetaStateElement,ThetaStateIterator,PhiElementParams
    EDCAInput = PiFactor,ThetaStateActionElement,TimeSteps,EDCAs,OtherAPs,ArrivalRate,InitialPacketInQueue
    # delay,StateActionMemo,SumNablaLogPi = zip(*[do_edca(EDCAInput) for i in range(Repeat)])
    delay,StateActionMemo,SumNablaLogPi = zip(*Parallel(n_jobs=-1)( [delayed(do_edca)(EDCAInput) for i in range(Repeat)] ))
    return np.sum(delay) / Repeat

def show_bar(Updates,Degrees,APs,OtherAPs):
    Proposed = [np.load('JThetas.npy')[Updates]]
    Conventional = [[[0,0,0],[3895,11495,20807],[12651,0,0]][APs][OtherAPs]]
    Heuristics = [[[0,0,0],[3804,10424,18298],[11507,0,0]][APs][OtherAPs]]
    Conventionaleft = [1]
    Heuristicleft = [2]
    Proposedleft = [3]

    width=0.8
    plt.bar(Proposedleft, Proposed, color='r', width=width, align='center')
    plt.bar(Conventionaleft, Conventional, color='g', width=width, align='center')
    plt.bar(Heuristicleft, Heuristics, color='b', width=width, align='center')

    center = [1,2,3]
    labels = ["Conventional","Heuristic","Proposed"]
    plt.xticks(center, labels)
    plt.ylim(10000,13000)

    plt.ylabel("Transmission delay time (μs)")
    plt.rcParams['ytick.direction'] = 'in'
    plt.savefig('bar.pdf')



def calculate_save_and_show(Updates,Episodes,TimeSteps,Degrees,EDCAs,APs,OtherAPs,LearningRate,ArrivalRate,InitialPacketInQueue,Repeat,mode,program):
    PhiElementParams = (np.ones(APs*5+OtherAPs*2,dtype=int),(1/5)*np.ones(APs*5+OtherAPs*2,dtype=int))
    # PhiElementParams = (np.ones(APs*5+OtherAPs*2,dtype=int),(1/5)*np.ones(APs*5+OtherAPs*2,dtype=int))
    # PhiElementParams[1][0:2] = 1/5
    for p in program:
        if p==0:
            print(delay_of_static_policy(Repeat,mode,EDCAs,APs,OtherAPs,TimeSteps,ArrivalRate,InitialPacketInQueue))
            # OtherAPs=0,mode=VO:    3895
            # OtherAPs=0,mode=less:  3804
            # OtherAPs=1,mode=VO:    11495
            # OtherAPs=1,mode=less:  10424
            # OtherAPs=2,mode=VO:    20807
            # OtherAPs=2,mode=less:  18298
        if p==1:
            calculate_and_save_result(Updates,Episodes,TimeSteps,Degrees,EDCAs,APs,OtherAPs,LearningRate,ArrivalRate,InitialPacketInQueue,PhiElementParams)
        if p==2:
            load_and_show_learning_curve(Updates,APs,OtherAPs)
        if p==3:
            load_and_show_policy(Updates,Degrees,EDCAs,APs,OtherAPs,PhiElementParams)
        if p==4:
            print(delay_of_gotten_policy(Repeat,Updates,Degrees,PhiElementParams,EDCAs,OtherAPs,TimeSteps,ArrivalRate,InitialPacketInQueue))
        if p==5:
            show_bar(Updates,Degrees,APs,OtherAPs)

calculate_save_and_show(Updates=100,Episodes=10000,TimeSteps=10,Degrees=3,EDCAs=2,APs=2,OtherAPs=0,LearningRate=1*10**-4,ArrivalRate=[1/200,1/400,0,0],InitialPacketInQueue=[0,0,0,0],Repeat=100000,mode='less',program=[3])

print(np.load('Jthetas.npy')[100])
print(100-100*np.load('Jthetas.npy')[0]/12651)
print(100-100*np.load('Jthetas.npy')[0]/11507)

# ////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////
# 送信遅延時間の計算
# ////////////////////////////////////////////////////////////////////////////////////////////////////
# ////////////////////////////////////////////////////////////////////////////////////////////////////
#
# Updates = 10
# Repeat = 10
# Episodes = 100
# Degrees = 3
# mode='VO'
# EDCAs=2
# APs = 2
# OtherAPs=0
# TimeSteps=10
# ArrivalRate=[1/200,1/300,1/400,1/500]
# InitialPacketInQueue=[0,0,0,0]
# LearningRate = 10**-4
# k=0
# PhiElementParams = (np.ones(OtherAPs*2+5,dtype=int),(1/5)*np.ones(OtherAPs*2+5,dtype=int))


if __name__ == '__main__':
    import time as ttiimmee
    start = ttiimmee.time()

    calculate_save_and_show(Updates=100,Episodes=10000,TimeSteps=10,Degrees=3,EDCAs=2,APs=3,OtherAPs=0,LearningRate=1*10**-4,ArrivalRate=[1/200,1/400,0,0],InitialPacketInQueue=[0,0,0,0],Repeat=1,mode='less',program=[1])

    elapsed_time = ttiimmee.time() - start
    print ("elapsed_time:{0}".format(elapsed_time) + "[sec]")

# elapsed_time:107.05730676651001[sec]
# elapsed_time:133.29407596588135[sec]

a = np.zeros(3,dtype=int)
a[0] = 2**60
a[0] += 2**60
print(a)




# joblibで更新
# npできるとこはnp
