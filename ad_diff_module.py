import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm , trange
import math
from qiskit import QuantumCircuit,execute,Aer,IBMQ,QuantumRegister,ClassicalRegister



##########################################その他関数################################################
def encoding_bin_num(qc,bin_num,bin_index):
    for i,_bin in enumerate(reversed(bin_num)):
        if _bin == '1':
            qc.x(bin_index[i])
            
def encoding_num(qc,num,index,qnum_after_decimal,qnum_interger):
    binnum = quantization(num,qnum_after_decimal,qnum_interger)
    for i,_bin in enumerate(reversed(binnum)):
        if _bin == '1':
            qc.x(index[i]) 
            
def encoding_const_2D(circuit,C,D,index_dic,qnum_after_decimal,qnum_interger):
    conD = quantization(D,qnum_after_decimal,qnum_interger)
    conDtC = quantization(D+C,qnum_after_decimal,qnum_interger)
    conDDDDCC = quantization(1-4*D-4*C,qnum_after_decimal,qnum_interger)
    encoding_bin_num(circuit,conD,index_dic['conD'])
    encoding_bin_num(circuit,conDtC,index_dic['conDtC'])
    encoding_bin_num(circuit,conDDDDCC,index_dic['conDDDDCC'])
    
    
def encoding_const_ippan(circuit,const,index,qnum_after_decimal,qnum_interger):
    conbin = quantization(const,qnum_after_decimal,qnum_interger)
    encoding_bin_num(circuit,conbin,index)
    
def encoding_const(circuit,C,D,index_dic,qnum_after_decimal,qnum_interger):
    bin_cons_x = quantization(1-2*D-C,qnum_after_decimal,qnum_interger)
    bin_cons_y = quantization(D,qnum_after_decimal,qnum_interger)
    bin_cons_z = quantization(D+C,qnum_after_decimal,qnum_interger)
    encoding_bin_num(circuit,bin_cons_x,index_dic['x_c'])
    encoding_bin_num(circuit,bin_cons_y,index_dic['y_c'])
    encoding_bin_num(circuit,bin_cons_z,index_dic['z_c'])
    
def encoding_three(circuit,qi_1,qi,qit1,index_dic,qnum_after_decimal,qnum_interger):
    bin_x = quantization(qi,qnum_after_decimal,qnum_interger)
    bin_y = quantization(qit1,qnum_after_decimal,qnum_interger)
    bin_z = quantization(qi_1,qnum_after_decimal,qnum_interger)
    encoding_bin_num(circuit,bin_x,index_dic['x'])
    encoding_bin_num(circuit,bin_y,index_dic['y'])
    encoding_bin_num(circuit,bin_z,index_dic['z'])


def execute_qc_mps_simulator(qc):   
    backend = provider.get_backend('simulator_mps')
    shots = 10
    results = execute(qc,backend=backend,shots=shots).result()
    answer = results.get_counts()
    return answer

def execute_qc(qc,provider,shots,simulator_name):   
    backend = provider.get_backend(simulator_name)
    results = execute(qc,backend=backend,shots=shots).result()
    answer = results.get_counts()
    return answer

# 量コン解を二進数に変換
def quant_ans_to_bin_ans(answer):
    for ans in answer:
        ans
    return ans


def ans_to_syousuu_sitei(binA,ika):
    k_jou = ika * (-1)
    sum = 0
    for i in reversed(binA):
        if i == '1':
            sum += 2**k_jou
        k_jou += 1
    return sum

# 量コンの解を整数へ変換する関数
def ans_to_int(answer):
    # 量コンからstr型binary取り出し
    for ans in answer:
        ans
    # 文字型二進数を数字型に進数に
    int_number = int(ans,2)
    return int_number

# #################### fourier / inv-fourier #####################
def qft(qc,num_index):
    for i in reversed(num_index):
        qc.h(i)
        k = 2
        for j in reversed(range(num_index[0],i)):
            qc.append(c_pk(k),[j,i])
            k += 1
        qc.barrier()
        
    
# 逆フーリエ変換
def qft_dagger(qc,num_index):
    control_k = 1
    for j in num_index:
        k = control_k
        for i in range(num_index[0],j):
            if i == j:
                continue
            qc.append(inv_c_pk(k),[i,j])
            k -= 1
        qc.h(j)
        qc.barrier()
        control_k += 1
    


# コントロールPゲート
def c_pk(k):
    P = QuantumCircuit(1)
    P.p(2 * math.pi/(2**k),0)
    p = P.to_gate()
    c_p = p.control(1)
    return c_p

# invコントロールPゲート
def inv_c_pk(k):
    P = QuantumCircuit(1)
    P.p(-2 * math.pi/(2**k),0)
    p = P.to_gate()
    c_p = p.control(1)
    return c_p



def qft_gate(_qnum):
    qft_circuit = QuantumCircuit(_qnum)
    qft(qft_circuit,[i for i in range(_qnum)])
    return qft_circuit.to_instruction()

def qft_dagger_gate(_qnum):
    qft_dagger_circuit = QuantumCircuit(_qnum)
    qft_dagger(qft_dagger_circuit,[i for i in range(_qnum)])
    return qft_dagger_circuit.to_instruction()
    

def fourier_swap(qc,index,qnum):
    for i in range(math.floor(qnum / 2)):
        qc.swap(index[i],index[-(i+1)])




################################### 量子化関数一覧
def quantization(num,qnum_after_decimal,qnum_interger):
    interger_part = math.floor(num)
    decimal_part = num - interger_part
    answer_decimal = ''
    differ = 1

    for i in [decimal_part_to_binary(i,qnum_after_decimal) for i in range(2**qnum_after_decimal)]:
        can_bin = ans_to_syousuu(i,qnum_after_decimal)
        sa = abs(decimal_part - can_bin)
        if differ > sa:
            answer_decimal = i
            differ = sa         
    answer_interger = interger_part_to_binary(interger_part,qnum_interger)
    return answer_interger + answer_decimal

def decimal_part_to_binary(float_num,qnum_after_decimal):
    sitei = '0{0}b'.format(qnum_after_decimal)
    bin_num = format(float_num,sitei)
    return bin_num

def interger_part_to_binary(int_num,qnum_interger):
    sitei = '0{0}b'.format(qnum_interger)
    bin_num = format(int_num,sitei)
    return bin_num

# バイナリーの解を決められた小数点以下のfloatに変換する
def ans_to_syousuu(binA,qnum_after_decimal):
    k_jou = qnum_after_decimal * (-1)
    sum = 0
    for i in reversed(binA):
        if i == '1':
            sum += 2**k_jou
        k_jou += 1
    return sum



###########################################足し算関数##########################################
# 加算部分
def addition(qc,num1_index,num2_index):    
    for j in reversed(range(0,len(num1_index))):
        k = 1
        for i in reversed(range(num1_index[0],num1_index[j]+1)):
            qc.append(c_pk(k),[i,num2_index[j]])
            k += 1
        qc.barrier()

# 足し算のゲート
def addition_gate(_qnum):
    add_circuit = QuantumCircuit(_qnum*2)
    addition(add_circuit,[i for i in range(_qnum)],[i for i in range(_qnum,_qnum*2)])
    return add_circuit.to_instruction()
    
# 足し単体
# 決められた数だけ足すゲート
# +1のときはtanni 1
def c_pk_tanni(tanni):
    add_x_circuit = QuantumCircuit(qnum)
    k = 0
    for i in range(qnum):
        add_x_circuit.p(2**(tanni-2) * 2 * math.pi/(2**k),i)
        k+=1
    add_x_circuit = add_x_circuit.to_gate()
    c_pk_tanni = add_x_circuit.control(1)
    return c_pk_tanni


# In[48]:


# ##########################################掛け算関数##########################################
# コントロールコントロールPゲート
def cc_pk(tanni,_qnum):
    P = QuantumCircuit(_qnum)
    for i in range(_qnum):
        P.p(2** i * tanni * 2 * math.pi/(2**_qnum),i)
    p = P.to_gate()
    cc_p = p.control(2)
    return cc_p

# 掛け算サーキット　インデックスの引数は、6桁9桁解入れる9桁
# swapを前後でつける必要あり。→内臓したい。

def multi_gate(_qnum,_qnum_const):
    multi_circuit = QuantumCircuit(_qnum*2 + _qnum_const)
    multi_circuit_list = [i for i in range(_qnum*2 + _qnum_const)]


    aindex = []
    bindex = []
    cindex = []
    for i in range(_qnum_const):
        aindex.append(multi_circuit_list[0])
        multi_circuit_list.pop(0)
    for i in range(_qnum):
        bindex.append(multi_circuit_list[0])
        multi_circuit_list.pop(0)
    cindex = multi_circuit_list


    control_tanni = 1
    for j in bindex:
        tanni = control_tanni
        for i in aindex:
            multi_circuit.append(cc_pk(tanni,_qnum),[i,j]+cindex)
            tanni *= 2
        control_tanni *= 2

    return multi_circuit.to_instruction()



# 1格子点の移流拡散計算をする回路 量子拡散回路

def ad_diff_circuit(qnum_dic,const_dic,index_dic):
    # 量子回路準備
    qnum = qnum_dic['qnum']
    qnum_const = qnum_dic['qnum_const']
    qnum_after_decimal = qnum_dic['qnum_after_decimal']
    qnum_interger = qnum_dic['qnum_interger']
    C = const_dic['C']
    D = const_dic['D']
    bin_x = QuantumRegister(qnum,'x')
    const_x = QuantumRegister(qnum_const,'con_x')
    bin_y = QuantumRegister(qnum,'y')
    const_y = QuantumRegister(qnum_const,'con_y')
    bin_z = QuantumRegister(qnum,'z')
    const_z = QuantumRegister(qnum_const,'con_z')
    answer = QuantumRegister(qnum ,'answer')
    qc = QuantumCircuit(bin_x,const_x,bin_y,const_y,bin_z,const_z,answer)

    encoding_const(qc,C,D,index_dic,qnum_after_decimal,qnum_interger) # 定数のエンコード
    qc.barrier()
    qc.append(qft_gate(qnum),index_dic['ans']) #フーリエ変換
    fourier_swap(qc,index_dic['ans'],qnum)
    qc.append(multi_gate(qnum,qnum_const),index_dic['x_c'] + index_dic['x'] + index_dic['ans'])  # (1-2D-C)×qi 
    qc.append(multi_gate(qnum,qnum_const),index_dic['y_c'] + index_dic['y'] + index_dic['ans']) #Dqi+1
    qc.append(multi_gate(qnum,qnum_const),index_dic['z_c'] + index_dic['z'] + index_dic['ans']) # (D+C)qi-1 
    # 逆フーリエ変換 & 測定
    fourier_swap(qc,index_dic['ans'],qnum)
    qc.append(qft_dagger_gate(qnum),index_dic['ans'])
    return qc.to_instruction()



#################################### 拡散1stepデータ受け取りからデータ出力まで ##############################
# 横軸xの格子点-2 回/step 計算が行われる


# def solve_ad_diff_1step(q_step,qnum,qnum_const,qnum_after_decimal,qnum_interger,index_dic,C,D,x,provider,shots,simulator):
#     measure_index = np.array([i for i in range(qnum)])
#     addiff_index = [i for i in range(qnum*4 + qnum_const * 3)]
#     circuit = QuantumCircuit(4*qnum + 4*qnum_const,qnum*(x-2))
#     add_diff_circuit = ad_diff_circuit(qnum,qnum_const,qnum_after_decimal,qnum_interger,C,D,index_dic)
    
#     for i in range(x-2):
#         encoding_three(circuit,q_step[i],q_step[i+1],q_step[i+2],index_dic,qnum_after_decimal,qnum_interger)
#         circuit.append(add_diff_circuit,addiff_index)
#         circuit.measure(index_dic['ans'],measure_index +i*qnum)
#         circuit.reset(addiff_index)
#         circuit.barrier()
#     answer = execute_qc(circuit,provider,shots,simulator)
    
#     answer_bin = quant_ans_to_bin_ans(answer)
#     ans = []
#     ans.insert(0,answer_bin[-qnum:])
#     index_st = -qnum
#     index_end = -qnum*2
#     for i in range(x-3):
#         ans.insert(0,answer_bin[index_end:index_st])
#         index_end-=(qnum)
#         index_st-=(qnum)

#     syousuu_ans = []
#     for i in ans:
#         syousuu_ans.insert(0,ans_to_syousuu_sitei(i,qnum_after_decimal*2))
        
#     # 0調整
#     syousuu_ans.append(0)
#     syousuu_ans.insert(0,0)
#     return syousuu_ans




def solve_addiff_1step(q_step,qnum_dic,index_dic,const_dic,execute_setting):
    ad_diff_circ = ad_diff_circuit(qnum_dic,const_dic,index_dic)

    for i in trange(10):
        wb = load_workbook(filename=file_path)
        ws = wb['Sheet1']
        data_list = []
        for j in excel_columns:
            data_list.append(ws[j+excel_index[i]].value)
        # 数合わせ
        data_list.append(0)
        data_list.insert(0,0)

        circ_list = []
        for k in range(const_dic['x']):
            length = qnum_dic['qnum']*4+qnum_dic['qnum_const']*3
            qc = QuantumCircuit(length,qnum_dic['qnum'])
            adm.encoding_three(qc,data_list[k],data_list[k+1],data_list[k+2],index_dic,qnum_dic['qnum_after_decimal'],qnum_dic['qnum_interger'])
            qc.append(ad_diff_circ,[_ for _ in range(length)])
            qc.measure(index_dic['ans'],[_ for _ in range(qnum_dic['qnum'])])
            circ_list.append(qc)

        # #リストで1ｓtep解が得られる。   
        quantum_answer = adm.execute_qc(circ_list, execute_setting['provider'], execute_setting['shots'], execute_setting['simulator'])
        new_data_list = [adm.ans_to_syousuu_sitei(adm.quant_ans_to_bin_ans(i),qnum_dic['qnum_after_decimal']*2) for i in quantum_answer]

        for j,column_num in enumerate(excel_columns):
            ws[column_num+excel_index[i+1]] = new_data_list[j]
        wb.save(file_path)