''' 
This file contains an implementation of the SIMPLER algorithm.

Written by:
    Rotem Ben-Hur - 
    Ronny Ronen - 
    Natan Peled - natanpeled@campus.technion.ac.il

Versions history:

The "RUN TIME parameters" part includes the next parameters:
    BenchmarkStrings - a list of minimized NOR and NOT netlists files to apply the algorithm on. The type of the files should be *.v.
    ROW_SIZE - a list of row sizes to run the algorithm with.
    JSON_CODE_GEN - to create an execution sequence JSON file, set this flag to TRUE.
    PRINT_CODE_GEN - to enable information print, set the flag to True. 
    PRINT_WARNING - to enable warnings print, set the flag to True.
    MAX_D - TODO - need to complete
    SORT_ROOTS - for arbitrary roots order set to 'NO'. For ascending roots order (by CU value) set to 'ASCEND'.
                 For descending roots order (by CU value) set to 'DESCEND'.

''' 

import numpy as np
import json


#====================== RUN TIME parameters =====================
#BenchmarkStrings = ['adder_netlist.v','arbiter_netlist.v','bar_netlist.v','cavlc_netlist.v','ctrl_netlist.v','dec_netlist.v','full_adder_12gates.v',
#                    'int2float_netlist.v','max_netlist.v','mult8-2n.v','priority_netlist.v','sin_netlist.v','voter_netlist.v'] #EPFL


#print controls
JSON_CODE_GEN = False
PRINT_CODE_GEN = False
PRINT_WARNING = True
SORT_ROOTS = 'NO' #Set to one of the follows: 'NO, 'ASCEND' 'DESCEND' 
 
#limit max D
#MAX_D = 20000 
#ROW_SIZE = [256,512,1024] #a list of row sizes to run the algorithm with


#================== End of RUN TIME parameters ==================


#================ Globals variables and Classes =================

class CellState:
    
    #States declarations 
    available = 1 #The cell was allocated and available again
    used = 2 #The cell is in use
    init = 3 #The cell is not in use, but need to initialized  
    
    #Global class parameters 
    max_num_of_used_cells = 0
    cur_num_of_used_cells = 0

    #Class methods:     
    @classmethod
    def Set_max_num_of_used_cells_to_zero(cls):
        cls.max_num_of_used_cells = 0
    
    @classmethod    
    def get_max_num_of_used_cells(cls):
        return cls.max_num_of_used_cells
     
    @classmethod    
    def Set_cur_num_of_used_cells_to_zero(cls):
        cls.cur_num_of_used_cells = 0
    
    @classmethod
    def get_cur_num_of_used_cells(cls):
        return cls.cur_num_of_used_cells

    @classmethod
    def State_available(cls):
        return cls.available
    
    @classmethod
    def State_used(cls):
        return cls.used
    
    @classmethod
    def State_init(cls):
        return cls.init 
       
    #End of class methods
    
    def __init__(self):
        self.state = CellState.available
        self.current_gate = -1
       
    def MarkAsAvailable(self):
        self.state = CellState.available

    def MarkAsUsed(self,gate_num):
        self.state = CellState.used
        self.current_gate = gate_num 
        CellState.cur_num_of_used_cells += 1
        if CellState.cur_num_of_used_cells > CellState.max_num_of_used_cells:
            CellState.max_num_of_used_cells = CellState.cur_num_of_used_cells
        
    def MarkAsInit(self):
        self.state = CellState.init
        CellState.cur_num_of_used_cells -= 1
        
    def GetCellState(self):
        return self.state
    
    def GetCurGate(self):
        return self.current_gate
        
    def Print(self):
        print('state is',self.state,'gate is',self.current_gate)
        
         
class Code_Generation_Table_Line:
    #Holds the data about a specific net, including allocated cells, allocations times and input's indexes list
    # In a case of INITIALIZATION, the line holds the initialized cells list instead of inputs list (other fields set to default values) 
    
    #Op declarations 
    No_inputs = 'NO INPUTS'
    Input = 'INPUT'
    Initialization = 'INITIALIZATION'
    
    #Class methods:        
    @classmethod
    def Get_no_inputs_op_val(cls):
        return cls.No_inputs
    
    @classmethod
    def Get_Input_op_val(cls):
        return cls.Inputs
    
    @classmethod
    def Get_Initialization_op_val(cls):
        return cls.Initialization
    
    #End of class methods
      
    def __init__(self,serial_number=0,inputs_list=[],op='',cell=0,time=0):
        self.serial_number = serial_number
        self.inputs_list = inputs_list
        self.op = op
        self.cell = cell
        self.time = time
    
    #Seters/geters:         
    def Set_serial_number(self,val):
        self.serial_number = val

    def Get_serial_number(self):
        return self.serial_number

    def Set_op(self,val):
        self.op = val

    def Get_op(self):
        return self.op
    
    def Set_cell(self,val):
        self.cell = val

    def Get_cell(self):
        return self.cell   
    
    def Set_time(self,val):
        self.time = val

    def Get_time(self):
        return self.time 
    
    def Insert_To_inputs_list(self,val):
        self.inputs_list.append(val)
        
    def Get_inputs_list(self):
        return self.inputs_list
    
    #End of seters/geters
    
    #Other methods:    
    def LineIntrlPrint(self):
        print('serial_number =',self.serial_number,'inputs_list =',self.inputs_list,'op =',self.op,'cell =',self.cell,'time =',self.time)
    
    def InsertInputLine(self,SN):
        #Updates the table - inserts a line who describes a net input
        self.serial_number = SN
        self.inputs_list = []
        self.op = self.Input        
        self.cell = SN
        self.time = 0
    
    def Insert_readoperations_parameters(self,SN,input_idxs,op):
        self.serial_number = SN
        self.inputs_list = input_idxs
        self.op = op

        
    def Insert_AllocateCell_parameters(self,cell,time):
        self.cell = cell
        self.time = time
        
    def Insert_No_Input_Line(self,SN):
        #Updates the table - inserts a line who describes a gate/wire who doesn't have inputs
        self.serial_number = SN
        self.inputs_list = []
        self.op = self.No_inputs        
        self.cell = 0
        self.time = 0 

# End of class Code_Generation_Table_Line   
           
class Code_Generation_Top_Data_struct:   
    #Holds all necessary data for statistics and JSON file creation
     
    def __init__(self,NetSizeRow,NetSizeCol,RowSize,NumberOfGates):
        self.code_generation_table = [Code_Generation_Table_Line() for idx in range(0,NetSizeRow)]
        self.TotalCycles = 0
        self.ReuseCycles = 0
        self.InitializationPercentage = 0.0
        self.NumberOfGates = NumberOfGates
        self.RowSize = RowSize
        self.InitializationList = [] #composed of Code_Generation_Table_Line instances 
        self.NoInputWireNum = 0
        self.NoInputWireList = []
        self.lr = NetSizeRow
        self.lc = NetSizeCol
        self.Max_Num_Of_Used_Cells = 0    
    
    #Seters/geters:       
    def Increase_ReuseCycles_by_one(self):
        self.ReuseCycles += 1 

    def Inset_To_InitializationList(self,val):
        self.InitializationList.append(val)
        
    def Get_InitializationList(self):
        return self.InitializationList
    
    def Increase_NoInputWireNum_by_one(self):
        self.NoInputWireNum += 1
        
    def Get_NoInputWireNum(self):
        return self.NoInputWireNum

    def Inset_To_NoInputWireList(self,val):
        self.NoInputWireList.append(val)
        
    def Get_NoInputWireList(self):
        return self.NoInputWireList   
    
    def Set_Max_Num_Of_Used_Cells(self,val):  
        self.Max_Num_Of_Used_Cells = val
    
    def Get_Max_Num_Of_Used_Cells(self):
        return self.Max_Num_Of_Used_Cells
    
    #End of Seters/geters 
    
    #Other methods:    
    #The next 4 methods are wrappers for the corresponding methods on Code_Generation_Table_Line class    
    def InsertInputLine(self,SN):
        self.code_generation_table[SN].InsertInputLine(SN)
        
    def Insert_readoperations_parameters(self,SN,input_idxs,op):
        self.code_generation_table[SN].Insert_readoperations_parameters(SN,input_idxs,op)
        
    def Insert_AllocateCell_parameters(self,line_num,cell,time):
        self.code_generation_table[line_num].Insert_AllocateCell_parameters(cell,time)
        
    def Insert_No_Input_Line(self,SN):
        self.code_generation_table[SN].Insert_No_Input_Line(SN) 
    #End of wrappers     
         
    def Add_To_Initialization_List(self,time,cells_list):
        #Adds an element to the Initialization
        
        self.Inset_To_InitializationList(Code_Generation_Table_Line(None,cells_list,Code_Generation_Table_Line.Get_Initialization_op_val(),None,time))
        self.Increase_ReuseCycles_by_one()          
    
    def Intrl_Print(self,Str):
        global PRINT_CODE_GEN
        if PRINT_CODE_GEN == True:
            print(Str)
                
    def PrintCodeGeneration(self,Benchmark_name):    
        #Prints the data and the statistics, can create the benchmark's execution sequence JSON file
        
        global varLegendRow,PRINT_CODE_GEN,InputString,OutputString,Benchmark
            
        #inputs
        input_list_for_print = 'Inputs:{' 
        for input_idx in range(0,self.lr - self.lc):
            if (self.code_generation_table[input_idx].Get_op() != Code_Generation_Table_Line.Get_no_inputs_op_val()): 
                input_list_for_print += InputString[input_idx] + '(' + str(self.code_generation_table[input_idx].Get_cell()) + '),'                    
        input_list_for_print = input_list_for_print[:len(input_list_for_print) - 1] + '}'
        self.Intrl_Print(input_list_for_print)
        
        #outputs
        output_len = len(OutputString)
        ofst = self.lr - output_len 
        output_list_for_print = 'Outputs:{' 
        for output_idx in range(0,output_len):
            idx = output_idx + ofst
            if (self.code_generation_table[idx].Get_op() != Code_Generation_Table_Line.Get_no_inputs_op_val()): #Has inputs
                output_list_for_print += OutputString[output_idx] + '(' + str(self.code_generation_table[idx].Get_cell()) + '),'                    
        output_list_for_print = output_list_for_print[:len(output_list_for_print) - 1] + '}'
        self.Intrl_Print(output_list_for_print)
                       
        #Execution sequence
        mergerd_table = self.code_generation_table + self.InitializationList
        mergerd_table.sort(key = lambda k: k.Get_time(), reverse=False) #Sorts by time
        execution_dict_for_JSON={} #JSON
        self.Intrl_Print('\nExecution sequence: {')
        for table_line in mergerd_table:
            if (table_line.Get_op() == Code_Generation_Table_Line.Get_Initialization_op_val()):            
                init_list_to_print = '{'
                for pair in table_line.Get_inputs_list(): #in a case of Initialization, inputs_list composed of [gate_number,cell_number] elements
                    init_list_to_print += varLegendRow[pair[0]] + '(' + str(pair[1]) + '),'
                init_list_to_print = init_list_to_print[:len(init_list_to_print) - 1] + '}' 
                self.Intrl_Print('T' + str(table_line.Get_time()) + ':Initialization(Ron)' +  init_list_to_print)
                execution_dict_for_JSON.update({'T' + str(table_line.Get_time()) : 'Initialization(Ron)' +  init_list_to_print})  #JSON 
            else:    
                table_line_name = varLegendRow[table_line.Get_serial_number()]
                if (table_line.Get_time() != 0): #not an input
                    inputs_str = ''        
                    for Input in table_line.Get_inputs_list():
                        inputs_str = inputs_str + varLegendRow[Input] + '(' + str(self.code_generation_table[Input].Get_cell()) + ')' + ','
                    inputs_str = '{' + inputs_str[:len(inputs_str) - 1] + '}'
                    self.Intrl_Print('T' + str(table_line.Get_time()) + ':' + table_line_name + '(' + str(table_line.Get_cell()) +')=' + table_line.Get_op() + inputs_str)
                    execution_dict_for_JSON.update({'T' + str(table_line.Get_time()) : table_line_name + '(' + str(table_line.Get_cell()) +')=' + table_line.Get_op() + inputs_str}) #JSON
                elif (table_line.Get_op() != Code_Generation_Table_Line.Get_no_inputs_op_val()): #line.time = 0 -- inputs
                    self.Intrl_Print('T' + str(table_line.Get_time()) + ':' + table_line_name + '(' + str(table_line.Get_serial_number()) +')=' + table_line.Get_op())             
        self.Intrl_Print('}')         
        
        #Statistics
        print '\nRESULTS AND STATISTICS:'
        print 'Benchmark:',Benchmark_name
        print 'Total number of cycles:',self.TotalCycles
        print 'Number of reuse cycles:',self.ReuseCycles
        self.InitializationPercentage = self.ReuseCycles/self.TotalCycles
        print 'Initialization Percentage:',self.InitializationPercentage
        connected_gates = self.NumberOfGates - self.NoInputWireNum
        print 'Number of gates:',connected_gates
        print 'Max number of used cells:',self.Max_Num_Of_Used_Cells
        print 'Row size (number of columns):',self.RowSize,'\n\n'
           
        #JSON creation
        if (JSON_CODE_GEN == True):
            top_JSON_dict = {'Benchmark':Benchmark,'Row size':self.RowSize,'Number of Gates':self.NumberOfGates,
                             'Inputs':input_list_for_print[len('Inputs:'):],'Number of Inputs':len(InputString),
                             'Outputs':output_list_for_print[len('Outputs:'):],'Number of Outputs':len(OutputString), 
                             'Total cycles':self.TotalCycles,'Reuse cycles':self.ReuseCycles,'Execution sequence' : execution_dict_for_JSON}                    
            with open('JSON_' + str(self.RowSize) + '_' + Benchmark[:Benchmark.find('.')] + '.JSON','w') as f:
                json.dump(top_JSON_dict,f,indent=4)
            f.close()                

              
    def PrintLines(self):
        for line in self.code_generation_table:
            line.LineIntrlPrint()
        
    def PrintLines_InitializationList(self):
        for line in self.InitializationList:
            line.LineIntrlPrint()


#End of class Code_Generation_Top_Data_struct            
            
#============ End of Globals variables and Classes ==============


#======================== Functions =============================

def readfield(tline,bmfId):
    #Parses the Inputs/Outputs/Wires declarations into a list
      
    FieldString = []
    while isinstance(tline,str):
        splited_tline = tline[tline.find(' '):].split(',')
        for i in list(filter(None,splited_tline)):
            FieldString.append(i.replace(' ','').replace(';','').replace('\n','').replace('\t',''))
        if (tline.find(';') != -1):
            break
        else:
            tline = bmfId.readline()
    return  list(filter(None,FieldString))


def readoperations(bmfId,varLegendRow,varLegendCol):
    #Parses the netlist assignments into a graph (matrix) 
           
    global D,lc,lr,code_generation_top,PRINT_WARNING 
    tline = bmfId.readline()
    while isinstance(tline,str):
        operands = []
        if tline.find('//') != -1: #Ignore comment lines at operands part
            tline = tline[0:tline.find('//')]
        if (tline.find("endmodule") != -1):
            return 
        elif ((tline.find("buf ") != -1) or (tline.find("zero ") != -1) or (tline.find("one ") != -1)):
            if (PRINT_WARNING == True):
                print("** Warning ** unsupported operation: \'" + tline.replace('\n','') + "\'\n")
        elif ((tline.find("inv") != -1) or (tline.find("nor") != -1)):
            
            idx_op = tline.find("nor")
            if (idx_op != -1):
                op = tline[idx_op:idx_op + len('nor') + 1] #op = "nor(#inputs)"
            else:
                op = 'inv1' 
            remain = tline
            num_of_op = 0
            while (remain.find('.') != -1):
                remain = remain[remain.find('.'):]
                remain = remain[remain.find('(') + 1:]
                open_bracket_idx = remain.find('(')
                close_bracket_idx = remain.find(')')
                #Gets the operands
                if (open_bracket_idx != -1) and (open_bracket_idx < close_bracket_idx):
                    second_close_bracket_idx = remain.find(')',close_bracket_idx + 1)
                    operands.append(remain[0:second_close_bracket_idx].replace(' ','').replace('\t',''))
                else:
                    operands.append(remain[0:close_bracket_idx].replace(' ','').replace('\t',''))
                num_of_op += 1
            out = operands[num_of_op - 1] #Gets the output operand
            outIdx = varLegendCol.index(out)
            input_idxs = []
            for k in range(0,num_of_op - 1):
                inIdx = varLegendRow.index(operands[k])
                D[inIdx][outIdx] = 1
                input_idxs.append(inIdx) #Gate inputs list
            code_generation_top.Insert_readoperations_parameters(outIdx + (lr -lc),input_idxs,op) #For statistics      
        tline = bmfId.readline()            


#----------------------- D illustration: ------------------------      
#              Out  0   1   2   3   4 . . . (n -number of inputs) 
#    In            w0   w1                        out2
#    0 - in0 
#    1 - in1
#    2
#    .
#    .
#    .
#    k - w0
#    .
#    .
#    .
#    n-1 - out1
#    n - out2


def GetRoots(): 
    #Returns the graph (D) roots. Also calculates the FO array.
    
    global D,FO,varLegendRow,UnConnected_wire,code_generation_top,len_input_and_wire
    roots = []
    for i in range(0,D.shape[0]):
        row_sum = D[i,:].sum()
        if (row_sum == 0):
            if (sum((ChildrenWithoutInputs(i))) != 0):
                roots.append(i)
            else:
                if (PRINT_WARNING == True):
                    print('** Warning **',varLegendRow[i],'has no input')
                code_generation_top.Increase_NoInputWireNum_by_one()
                code_generation_top.Inset_To_NoInputWireList(i)
                code_generation_top.Insert_No_Input_Line(i)
        else:
            FO[i] = row_sum   
    print('\n')
    return roots
    
def GetParents(V_i):
    global lr,lc
    return (np.where(D[V_i,:]==1)[0] + (lr - lc)) #(lr - lc) is the vertex index offset between columns and rows (relevant only for wires and outputs)

def GetChildrens(V_i):
    return np.where(D[:,V_i]==1)[0]


def ChildrenWithoutInputs(V_i):
    #Returns the childrens without netlist inputs
    
    global LEAFS_inputs,lr,lc
    childrens = GetChildrens(V_i - (lr - lc)) #(lr - lc) is the vertex index offset between columns and rows (relevant only for wires and outputs)
    childrens_without_inputs = [child for child in childrens if (child in LEAFS_inputs) == False]
    return childrens_without_inputs
    
        
def computeCU(V_i):
    #Computes the cell usage (CU) of gate (node) V_i  
    
    global CU
    if (CU[V_i] > 0):
        return # CU[V_i] was already generated and therefore doesn't change   
    childrens = ChildrenWithoutInputs(V_i)          
    if(len(childrens) == 0): #V_i has no childrens -> V_i is connected to function inputs only
        CU[V_i] = 1
    else:
        if (len(childrens) == 1):
            computeCU(childrens[0])
            CU[V_i] = CU[childrens[0]]
        else:
            childrens_cu = []
            for child in childrens:
                computeCU(child)
                childrens_cu.append(CU[child])
            childrens_cu.sort(key=None, reverse=True)               
            num_of_childrens = range(0,len(childrens)) # equal to + (i - 1) for all i in 1 to N (N is the number of childrens)           
            CU[V_i] = max(np.add(childrens_cu,num_of_childrens))


def AllocateRow(V_i):
    #Allocates cells to the gate V_i and his children (a sub-tree rooted by V_i). 
    #In a case the allocation for one of V_i's children or V_i itself is failed, the function returns False. On successful allocation returns True.
    
    global Map
    childrens = ChildrenWithoutInputs(V_i) #Equal to C(V_i) - the set of V_i's childrens 
    childrens_sorted_by_cu = [[child,CU[child]] for child in childrens] # the loop creates a list composed of pairs of the form [child number, CU[child number]]
    childrens_sorted_by_cu.sort(key = lambda k: k[1], reverse=True) #sorting by CU   
    childrens_sorted_by_cu = [elm[0] for elm in childrens_sorted_by_cu] #taking only the child (vertex) number       
    for V_j in childrens_sorted_by_cu:
        if (Map[V_j] == 0):
            if (AllocateRow(V_j) == False):
                return False
    if (Map[V_i] == 0): #V_i is not mapped
        Map[V_i] = AllocateCell(V_i)         
        if (Map[V_i] == 0): #V_i could not be mapped
            return False
    return True


def AllocateCell(V_i):
    #Allocates cell for 1 gate. In case available cells don't exist, the function will initialize all the cells whose state is init, and will look again for available cell.
    #At this point, if an available cell still doesn't exist, the function will return 0 (there is no mapping).  
    
    global CELLS,lr,lc,i,N,t,FO,Map,code_generation_top
    CellsForInit = []
    FreeCell = 0
    j = 0 #To make the variable scope global in this function
    for j in range(i,N): #The N-i dedicated cells (columns) for the computation. j gets values in the range i<=j<=N
        #if j < N: #Only N cells exist  
        if (CELLS[j].GetCellState() == CellState.State_available()):
            FreeCell = j 
            break #Available cell exists 
        elif (CELLS[j].GetCellState() == CellState.State_init()):
            CellsForInit.append([CELLS[j].GetCurGate(),j]) # CELLS[j].current_gate is for Json file creation
    if (FreeCell == 0): #No available cell; therefore, all init cells are initialized simultaneously. 
        if (CellsForInit == []):
            return 0 #No cells to initialize 
        else:
            for j in range(N - 1,i - 1,-1): # at the for loop end's j is the first cell in CellsForInit 
                if (CELLS[j].GetCellState() == CellState.State_init()):
                    CELLS[j].MarkAsAvailable()
                    FreeCell = j
            t += 1 #Increments the number of cycles 
            code_generation_top.Add_To_Initialization_List(t,CellsForInit)            
    CELLS[FreeCell].MarkAsUsed(V_i) #Use the available cell
    t += 1 #Increments the number of cycles
    code_generation_top.Insert_AllocateCell_parameters(V_i,FreeCell,t)
    for V_k in ChildrenWithoutInputs(V_i): #Updates the FO value if the allocated cells
        FO[V_k] -= 1
        if (FO[V_k] == 0):
            CELLS[Map[V_k]].MarkAsInit()
    return FreeCell       
        


def IncreaseOutputsFo():
    # This function created to make sure outputs cells will not evacuated.
    #It is done by increasing their FO by 1
    
    global FO,len_input_and_wire,lr
    for idx in range(len_input_and_wire,lr): #outputs idx range
        FO[idx] += 1

#======================= End of Functions =======================    

   
#======================== SIMPLER MAPPING =======================
def SIMPLER_Main (BenchmarkStrings, MAX_D, ROW_SIZE, Benchmark_name):
    
    global CELLS,lr,lc,i,N,t,FO,Map,code_generation_top,len_input_and_wire,LEAFS_inputs,CU,D,PRINT_WARNING,varLegendRow,UnConnected_wire,varLegendRow,PRINT_CODE_GEN,InputString,OutputString,Benchmark,PRINT_CODE_GEN
    
    for Row_size in ROW_SIZE: 
        for Benchmark in BenchmarkStrings:
            
            #Parse operations 
            print("\n\n" + Benchmark + ". Started Benchmark " + str((BenchmarkStrings.index(Benchmark)) + 1) + "\n")
            bmfId = open(Benchmark,"r") #open file       
            #read input/output/wire
            tline = bmfId.readline()
            while isinstance(tline,str):
                comment_idx = tline.find('//')
                if (comment_idx != -1):
                    tline = tline[0:comment_idx - 1]
                input_idx = tline.find('input');
                output_idx = tline.find('output');
                wire_idx = tline.find('wire');
                if (input_idx != -1):
                    InputString = readfield(tline[input_idx:],bmfId)
                elif (output_idx != -1):
                    OutputString = readfield(tline[output_idx:],bmfId)
                elif (wire_idx != -1):
                    WireString = readfield(tline[wire_idx:],bmfId)
                    break
                tline = bmfId.readline()
            
            #Map variables to numbers
            varLegendCol = WireString + OutputString
            varLegendRow = InputString + WireString + OutputString 
            len_input_and_wire = len(InputString + WireString)
          
            gate_num = len(varLegendCol)
        
            #analyze dependencies
            lr = len(varLegendRow)
            lc = len(varLegendCol) 
            if (lr>MAX_D or lc>MAX_D):
                print("** net too big, skip " + str(lr) +" X " + str(lc) + "\n")
                continue
            code_generation_top = Code_Generation_Top_Data_struct(lr,lc,Row_size,gate_num) #creates a Code_Generation_Top_Data_struct instance       
            for input_idx in range(0,len(InputString)):   
                code_generation_top.InsertInputLine(input_idx)
                                      
            D = np.zeros((lr,lc),dtype = np.int) #Graph matrix
            readoperations(bmfId,varLegendRow,varLegendCol) #parses the netlist
            LEAFS_inputs = list(range(0,lr - lc)) #Inputs indexes
            UnConnected_wire = 0 #Number of unconnected wires (gates)
            code_generation_success_flag = True 
            
            #================ SIMPLER algorithm Starts ================ 
            
            FO = np.zeros(lr,dtype = np.int) #FO array initialization       
            ROOTs = GetRoots() #set of all roots of the graph. Also calculates FO values.
            IncreaseOutputsFo() # To ensure outputs who are also inputs, will not be evacuated
            t = 0 #Number of clock cycles
            
            CU = np.zeros(lr,dtype = np.int) #CU(V_i) is the Cell Usage of V_i. Array initialization 
            Map = np.zeros(lr,dtype = np.int) #Map(V_i) is the number of the cell/column V_i is mapped to. Array initialization
            
            i = lr - lc  #number of inputs
            N = Row_size
            CELLS = [CellState() for cell in range(0,N)] #CELLS list initialization
            for j in range(0,N):   
                if (j < i): # j = 0 -> (i - 1)
                    CELLS[j].MarkAsUsed(j) #an input
                else: # j = i -> N (N is equal to ROW_SISE)
                    CELLS[j].MarkAsAvailable()       
          
            
            for r in ROOTs:
                computeCU(r)      
            if SORT_ROOTS == 'NO':
                for r in ROOTs:    
                    if (AllocateRow(r) == False):
                        print('False - no mapping')
                        code_generation_success_flag = False #Printing flag 
                        break #To enable multiple runs. To fit the code to the article, comment this line, and uncomment the two next lines. 
                        #return False #Cannot find mapping 
                #return True #A mapping of the entire netlist was found 
            else:        
                sorted_ROOTs = [[r,CU[r]] for r in ROOTs]
                if SORT_ROOTS == 'DESCEND':
                    sorted_ROOTs.sort(key = lambda k: k[1], reverse=True)
                elif SORT_ROOTS == 'ASCEND':
                    sorted_ROOTs.sort(key = lambda k: k[1], reverse=False)
                for sr in sorted_ROOTs:    
                    if (AllocateRow(sr[0]) == False):
                        print('False - no mapping')
                        code_generation_success_flag = False
                        break
                        #return False #cannot find mapping
                #return True        
            
            #================ SIMPLER algorithm ends ================ 
                   
            #Statistics calculations 
            code_generation_top.TotalCycles = t
            code_generation_top.Set_Max_Num_Of_Used_Cells(CellState.get_max_num_of_used_cells())
            if (code_generation_success_flag == True):
                code_generation_top.PrintCodeGeneration(Benchmark_name) 
            
            #Benchmark's end 
            bmfId.close() #close file
            CellState.Set_cur_num_of_used_cells_to_zero() #need to initiate because its a class variable
            CellState.Set_max_num_of_used_cells_to_zero() #need to initiate because its a class variable
            print('End of Benchmark '+ str((BenchmarkStrings.index(Benchmark)) + 1) + '\n\n')
            
#=========================== End of SIMPLER MAPPING ===========================

#============================== End of code ===================================