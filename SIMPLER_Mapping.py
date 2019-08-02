''' 
This file contains an implementation of the SIMPLER algorithm.

Written by:
    Rotem Ben-Hur - rotembenhur@campus.technion.ac.il
    Ronny Ronen - 
    Natan Peled - natanpeled@campus.technion.ac.il

The "RUN TIME parameters" part includes the next parameters:
    ROW_SIZE - a list of row sizes to run the algorithm with.
    JSON_CODE_GEN - to create an execution sequence JSON file, set this flag to TRUE.
    PRINT_CODE_GEN - to enable information print, set the flag to True. 
    PRINT_WARNING - to enable warnings print, set the flag to True.
    Max_num_gates - the maximum number of gates the tool generates a mapping to 
    SORT_ROOTS - for arbitrary roots order set to 'NO'. For ascending roots order (by CU value) set to 'ASCEND'.
                 For descending roots order (by CU value) set to 'DESCEND'.

''' 

import numpy as np
import simplejson
from collections import OrderedDict
import time

#================ Globals variables and Classes =================


class GraphEdge:

    def __init__(self,Source,Dest,Val):
        self.s = Source
        self.d = Dest
        self.v = Val

    def GetDest(self):
        return self.d

    def GetSource(self):
        return self.s

    def GetVal(self):
        return self.v

    def Print(self):
        print('(' + str(self.s) + ',' + str(self.d) + '), edge val is:' + str(self.v))

class NodeData:
    
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
    
    def __init__(self,Num,Op='',Inputs_list=[],Time=0):
        self.node_num = Num 
        self.op = Op
        self.val = False
        self.FO = 0
        self.CU = 0
        self.map = 0
        self.inputs_list = Inputs_list
        self.time = Time
        self.SIMPLER_lists_node = None

        self.out_edges = []
        self.in_edges = []
        self.num_of_out_edges = 0
        self.num_of_in_edges = 0


    def SetNodeNum(self,Num):
        self.node_num = Num
    
    def GetNodeNum(self):
        return self.node_num
        
    def SetNodeOp(self,Op):
        self.op = Op
        
    def GetNodeOp(self):
        return self.op
        
    def SetNodeCu(self,Val):
        self.CU = Val
    
    def GetNodeCu(self):
        return self.CU

    def SetNodeFO(self,Val):
        self.FO = Val
    
    def GetNodeFO(self):
        return self.FO
        
    def SetNodeMap(self,Val):
        self.map = Val
        
    def GetNodeMap(self):
        return self.map
    
    def GetNodeTime(self):
        return self.time

    def SetNodeInputs_list(self,Inputs):
        self.inputs_list = Inputs
        
    def GetNodeInputs_list(self):
        return self.inputs_list        

    def InsertInputNode(self):
        #Creates a node who describes a net input
        self.inputs_list = []
        self.op = self.Input        
        self.map = self.node_num
        self.time = 0

    def Insert_readoperations_parameters(self,NodeNum,InputIdxs,Op):
        self.node_num = NodeNum
        self.inputs_list = InputIdxs
        self.op = Op
        
    def Insert_AllocateCell_parameters(self,cell,time):
        self.time = time
        self.map = cell
        
    def Insert_No_Input_Node(self,NodeNum):
        #Inserts a line who describes a gate/wire who doesn't have inputs
        self.node_num = NodeNum
        self.inputs_list = []
        self.op = self.No_inputs        
        self.map = -1 
        self.time = 0   

    def Set_SIMPLER_lists_node(self,ref):
        print(ref)
        self.SIMPLER_lists_node = ref

    def Get_SIMPLER_lists_node(self):
        return self.SIMPLER_lists_node
        
    def PrintNodeData(self):
        print('node_num =',self.node_num,'inputs_list =',self.inputs_list,'op =',self.op,'cell =',self.map,'time =',self.time,'CU =',self.CU,'FO =',self.FO)

    def AddOutEdge(self,Dest,Val):
        self.out_edges.append(GraphEdge(self.node_num,Dest,Val))
        self.num_of_out_edges += 1

    def AddInEdge(self,Source,Val):
        self.in_edges.append(GraphEdge(Source,self.node_num,Val))
        self.num_of_in_edges += 1
        
    def RemoveOutEdge(self,Dest):
        for i,e in enumerate(self.out_edges):
            if e.GetDest() == Dest:
                del self.out_edges[i]
                self.num_of_out_edges -= 1
                break

    def RemoveInEdge(self,Source):
        for i,e in enumerate(self.in_edges):
            if e.GetSource() == Source:
                del self.in_edges[i]
                self.num_of_in_edges -= 1
                break

    def GetNumOfOutEdges(self):
        return self.num_of_out_edges

    def GetNumOfInEdges(self):
        return self.num_of_in_edges

    def GetOutEdgesList(self,TrueDepFlag):
        if TrueDepFlag:
            return [e.GetDest() for e in self.out_edges if e.GetVal()==1]
        else:
            return [e.GetDest() for e in self.out_edges]

    def GetInEdgesList(self,TrueDepFlag):
        if TrueDepFlag:
            L = [e.GetSource() for e in self.in_edges if e.GetVal()==1]
            L.reverse()
            return L
        else:
            return [e.GetSource() for e in self.in_edges]

class CellInfo:
    
    #States declarations 
    Available = 1 #The cell was allocated and available again
    Used = 2 #The cell is in use
    Init = 3 #The cell is not in use, but need to initialized  
    
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

    #End of class methods
    
    def __init__(self,Idx):
        self.state = CellInfo.Available
        #self.current_gate = -1
        self.next = None
        self.prev = None
        self.cell_idx = Idx
        self.current_gate_num = None

    def GetCellIdx(self):
        return self.cell_idx

    def SetNext(self,idx):
        self.next = idx
        
    def SetPrev(self,idx):
        self.prev = idx

    def GetNext(self):
        return self.next
        
    def GetPrev(self):
        return self.prev

    def SetCurGateNum(self,gate_num):
        self.current_gate_num = gate_num

    def GetCurGateNum(self):
        return self.current_gate_num



class CellsInfo:

    def __init__(self,N):
        self.used_head = None
        self.used_tail = None
        self.init_head = None
        self.init_tail = None
        self.available_head = None
        self.available_tail = None
        self.cells = [CellInfo(idx) for idx in range(0,N)]
        self.init_list_for_json = []

    #Available list methods
    def GetFirst_Available(self):
        if self.available_head == None:
            return None
        else:
            return self.cells[self.available_head].GetCellIdx()

    def Concatenate_init_to_available_list(self):
        self.available_head = self.init_head
        self.available_tail = self.init_tail
        self.cells[self.available_tail].SetNext(None)
        self.cells[self.available_head].SetPrev(None)
        self.init_head = None
        self.init_tail = None

    def DeleteFirst_Available(self):
        if self.available_head == self.available_tail:
            self.cells[self.available_head].SetNext(None)
            self.cells[self.available_head].SetPrev(None)
            self.available_head = self.available_tail = None
        else:
            next_available_cell = self.cells[self.available_head].GetNext()
            self.cells[self.available_head].SetNext(None)
            self.available_head = next_available_cell
            self.cells[next_available_cell].SetPrev(None)

    def Insert_Available(self,cell_idx):
        if self.available_head == None:
            self.available_head = cell_idx
            self.available_tail = cell_idx
        else:
            self.cells[self.available_head].SetPrev(cell_idx)
            self.cells[cell_idx].SetNext(self.available_head)
            self.available_head = cell_idx
    #Init list methods
    def IsNotEmpty_Init(self):
        if self.init_head == None:
            return False
        else:
            return True

    def Empty_Init(self):
        self.init_tail = None
        self.init_head = None
        self.init_list_for_json = []

    def Insert_Init(self,cell_idx):
        if self.init_head == None:
            self.init_head = cell_idx
            self.init_tail = cell_idx
        else:
            self.cells[self.init_head].SetPrev(cell_idx)
            self.cells[cell_idx].SetNext(self.init_head)
            self.init_head = cell_idx
        cur_gate = self.cells[cell_idx].GetCurGateNum()
        self.init_list_for_json.append([cur_gate,cell_idx])

    #Used list methods
    def Insert_Used(self,cell_idx,gate_num):
        if self.used_head == None:
            self.used_head = cell_idx
            self.used_tail = cell_idx
        else:
            self.cells[self.used_head].SetPrev(cell_idx)
            self.cells[cell_idx].SetNext(self.used_head)
            self.used_head = cell_idx
        self.cells[cell_idx].SetCurGateNum(gate_num)

    def Delete_Used(self,cell_idx):
        next_cell =  self.cells[cell_idx].GetNext()
        prev_cell =  self.cells[cell_idx].GetPrev()

        if self.used_head == self.used_tail and self.used_tail == cell_idx:
            self.cells[cell_idx].SetNext(None)
            self.cells[cell_idx].SetPrev(None)
            self.used_tail = self.used_head = None
        elif cell_idx == self.used_head:
            next_head = self.cells[cell_idx].GetNext()
            self.cells[cell_idx].SetNext(None)
            self.used_head = next_head
        elif cell_idx == self.used_tail:
            next_tail = self.cells[cell_idx].GetPrev()
            self.cells[cell_idx].SetPrev(None)
            self.used_tail = None
        else:
            next_cell = self.cells[cell_idx].GetNext()
            prev_cell = self.cells[cell_idx].GetPrev()
            self.cells[cell_idx].SetNext(None)
            self.cells[cell_idx].SetPrev(None)
            self.cells[next_cell].SetPrev(prev_cell)
            self.cells[prev_cell].SetPrev(next_cell)
        
# End of class CellState 



class GraphNode:

    def __init__(self,NodeNum):
        self.node_num = NodeNum
        self.out_edges = []
        self.in_edges = []
        self.num_of_out_edges = 0
        self.num_of_in_edges = 0

    def AddOutEdge(self,Dest,Val):
        self.out_edges.append(GraphEdge(self.node_num,Dest,Val))
        self.num_of_out_edges += 1

    def AddInEdge(self,Source,Val):
        self.in_edges.append(GraphEdge(Source,self.node_num,Val))
        self.num_of_in_edges += 1
        
    def RemoveOutEdge(self,Dest):
        for i,e in enumerate(self.out_edges):
            if e.GetDest() == Dest:
                del self.out_edges[i]
                self.num_of_out_edges -= 1
                break

    def RemoveInEdge(self,Source):
        for i,e in enumerate(self.in_edges):
            if e.GetSource() == Source:
                del self.in_edges[i]
                self.num_of_in_edges -= 1
                break

    def GetNumOfOutEdges(self):
        return self.num_of_out_edges

    def GetNumOfInEdges(self):
        return self.num_of_in_edges

    def GetOutEdgesList(self,TrueDepFlag):
        if TrueDepFlag:
            return [e.GetDest() for e in self.out_edges if e.GetVal()==1]
        else:
            return [e.GetDest() for e in self.out_edges]

    def GetInEdgesList(self,TrueDepFlag):
        if TrueDepFlag:
            L = [e.GetSource() for e in self.in_edges if e.GetVal()==1]
            L.reverse()
            return L
        else:
            return [e.GetSource() for e in self.in_edges]

class SIMPLER_Top_Data_Structure:

    def __init__(self, RowSize, bmfId, Benchmark):
        self.PRINT_WARNING = False
        self.PRINT_CODE_GEN = False
        self.JSON_CODE_GEN = False    
        self.Benchmark = Benchmark
        self.InputString = []
        self.OutputString = []
        self.WireString = []
        self.varLegendRow = []
        self.varLegendCol = []
        self.len_input_and_wire = 0
        self.lr = 0
        self.lc = 0   
        self.RowSize = RowSize
        self.NumberOfGates = 0
        self.NodesList = []
        self.i = 0
        self.N = RowSize
        self.t = 0 #TotalCycles
        self.ReuseCycles = 0
        self.LEAFS_inputs = []
        self.NodesList = []
        self.InitializationList = [] #composed of NoedData dummy instances
        self.InitializationPercentage = 0.0
        self.NoInputWireNum = 0
        self.NoInputWireList = [] 
        self.Max_Num_Of_Used_Cells = 0  
        self.UnConnected_wire = 0
        self.cells = CellsInfo(self.N)


        # read input/output/wire
        tline = bmfId.readline()
        while isinstance(tline, str):
            comment_idx = tline.find('//')
            if (comment_idx != -1):
                tline = tline[0:comment_idx - 1]
            input_idx = tline.find('input');
            output_idx = tline.find('output');
            wire_idx = tline.find('wire');
            if (input_idx != -1):
                self.InputString = self.readfield(tline[input_idx:], bmfId)
            elif (output_idx != -1):
                self.OutputString = self.readfield(tline[output_idx:], bmfId)
            elif (wire_idx != -1):
                self.WireString = self.readfield(tline[wire_idx:], bmfId)
                break
            tline = bmfId.readline()
        # Map variables to numbers
        self.varLegendCol = self.WireString + self.OutputString
        self.varLegendRow = self.InputString + self.WireString + self.OutputString
        self.len_input_and_wire = len(self.InputString + self.WireString)
        self.lr = len(self.varLegendRow)
        self.lc = len(self.varLegendCol)   
        self.i = self.lr - self.lc  # number of inputs
        self.NodesList = [NodeData(idx) for idx in range(self.lr)]
        for input_idx in range(len(self.InputString)):   
            self.NodesList[input_idx].InsertInputNode()
        self.readoperations(bmfId)  # parses the netlist         
        self.LEAFS_inputs = list(range(self.i)) #Inputs indexes
        
    #Seters/geters:     
    def Get_lr(self):
        return self.lr
    
    def Get_lc(self):
        return self.lc
    
    def GetTotalCycles(self):
        return self.t
    
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

    #The next 4 methods are wrappers for the corresponding methods on Code_Generation_Table_Line class    
    def InsertInputNode(self,SN):
        self.NodesList[SN].InsertInputNode(SN)
        
    def Insert_readoperations_parameters(self,SN,input_idxs,op):
        self.NodesList[SN].Insert_readoperations_parameters(SN,input_idxs,op)
        
    def Insert_AllocateCell_parameters(self,line_num,cell):
        self.NodesList[line_num].Insert_AllocateCell_parameters(cell,self.t)
        
    def Insert_No_Input_Node(self,SN):
        self.NodesList[SN].Insert_No_Input_Node(SN) 
    #End of wrappers     

    def Add_To_Initialization_List(self,time,cells_list):
        #Adds an element to the Initialization
        self.Inset_To_InitializationList(NodeData(None,NodeData.Get_Initialization_op_val(),cells_list,time))
        self.Increase_ReuseCycles_by_one()          
            
    def Intrl_Print(self,Str):
        global PRINT_CODE_GEN
        if PRINT_CODE_GEN == True:
            print(Str)

    def PrintCodeGeneration(self):    
        #Prints the data and the statistics, can create the benchmark's execution sequence JSON file
                
        print('\\\\\\\\\\\\ MAPPING OF',self.Benchmark,'WITH ROW SIZE =',self.RowSize,' \\\\\\\\\\\\\n')
        
        #inputs
        input_list_for_print = 'Inputs:{' 
        for input_idx in range(0,self.lr - self.lc):
            if (self.NodesList[input_idx].GetNodeOp() != NodeData.Get_no_inputs_op_val()): 
                input_list_for_print += self.InputString[input_idx] + '(' + str(self.NodesList[input_idx].GetNodeMap()) + '),'                    
        input_list_for_print = input_list_for_print[:len(input_list_for_print) - 1] + '}'
        self.Intrl_Print(input_list_for_print)
        
        #outputs
        output_len = len(self.OutputString)
        ofst = self.lr - output_len 
        output_list_for_print = 'Outputs:{' 
        for output_idx in range(0,output_len):
            idx = output_idx + ofst
            if (self.NodesList[idx].GetNodeOp() != NodeData.Get_no_inputs_op_val()): #Has inputs
                output_list_for_print += self.OutputString[output_idx] + '(' + str(self.NodesList[idx].GetNodeMap()) + '),'                    
        output_list_for_print = output_list_for_print[:len(output_list_for_print) - 1] + '}'
        self.Intrl_Print(output_list_for_print)

        #Execution sequence
        mergerd_list = self.NodesList + self.InitializationList
        mergerd_list.sort(key = lambda k: k.GetNodeTime(), reverse=False) #Sorts by time
        #execution_dict_for_JSON=OrderedDict({}) #JSON
        cells_to_init_list = ['INIT_CYCLE(' + str(idx) + ')' for idx in range(self.lr-self.lc,self.lr)]
        execution_dict_for_JSON=OrderedDict({'T0':'Initialization(Ron)'+ str(cells_to_init_list).replace('[','{').replace(']','}').replace(' ','')}) #JSON
        self.Intrl_Print('\nEXECUTION SEQUENCE + MAPPING: {')
        for node in mergerd_list:
            if (node.GetNodeOp() == NodeData.Get_Initialization_op_val()):            
                init_list_to_print = '{'
                for pair in node.GetNodeInputs_list(): #in a case of Initialization, inputs_list composed of [gate_number,cell_number] elements
                    init_list_to_print += self.varLegendRow[pair[0]] + '(' + str(pair[1]) + '),'
                init_list_to_print = init_list_to_print[:len(init_list_to_print) - 1] + '}' 
                self.Intrl_Print('T' + str(node.GetNodeTime()) + ':Initialization(Ron)' +  init_list_to_print)
                execution_dict_for_JSON.update({'T' + str(node.GetNodeTime()) : 'Initialization(Ron)' +  init_list_to_print})  #JSON 
            else:    
                node_name = self.varLegendRow[node.GetNodeNum()]
                if (node.GetNodeTime() != 0): #not an input
                    inputs_str = ''        
                    for Input in node.GetNodeInputs_list():
                        inputs_str = inputs_str + self.varLegendRow[Input] + '(' + str(self.NodesList[Input].GetNodeMap()) + ')' + ','
                    inputs_str = '{' + inputs_str[:len(inputs_str) - 1] + '}'
                    self.Intrl_Print('T' + str(node.GetNodeTime()) + ':' + node_name + '(' + str(node.GetNodeMap()) +')=' + node.GetNodeOp() + inputs_str)
                    execution_dict_for_JSON.update({'T' + str(node.GetNodeTime()) : node_name + '(' + str(node.GetNodeMap()) +')=' + node.GetNodeOp() + inputs_str}) #JSON
                elif (node.GetNodeOp() != NodeData.Get_no_inputs_op_val()): #line.time = 0 -- inputs
                    self.Intrl_Print('T' + str(node.GetNodeTime()) + ':' + node_name + '(' + str(node.GetNodeNum()) +')=' + node.GetNodeOp())    
                    #execution_dict_for_JSON.update({'T' + str(node.GetNodeTime()) : node_name + '(' + str(node.GetNodeMap()) +')=' + node.GetNodeOp()}) #JSON                    
        self.Intrl_Print('}')         
        
        #Statistics
        print ('\nRESULTS AND STATISTICS:')
        print ('Benchmark:',self.Benchmark)
        print ('Total number of cycles:',self.t)
        print ('Number of reuse cycles:',self.ReuseCycles)
        self.InitializationPercentage = self.ReuseCycles/self.t
        print ('Initialization percentage:',self.InitializationPercentage)
        connected_gates = self.lc - self.NoInputWireNum
        print ('Number of gates:',connected_gates)
        #print ('Max number of used cells:',self.Max_Num_Of_Used_Cells)
        print ('Row size (number of columns):',self.RowSize,'\n\n')

        #JSON creation
        if (JSON_CODE_GEN == True):
            top_JSON_dict=OrderedDict({'Benchmark':self.Benchmark}) #JSON
            top_JSON_dict.update({'Row size':self.RowSize})
            top_JSON_dict.update({'Number of Gates':connected_gates})
            top_JSON_dict.update({'Inputs':input_list_for_print[len('Inputs:'):]})
            top_JSON_dict.update({'Outputs':output_list_for_print[len('Outputs:'):]})
            top_JSON_dict.update({'Number of Inputs':len(self.InputString)})
            top_JSON_dict.update({'Total cycles':self.t})
            top_JSON_dict.update({'Reuse cycles':self.ReuseCycles})
            top_JSON_dict.update({'Execution sequence' : execution_dict_for_JSON})
            #print('Benchmark= ',Benchmark)                             
            with open('JSON_' + str(self.RowSize) + '_' + self.Benchmark + '.json','w') as f:
                simplejson.dump(top_JSON_dict,f,indent=4)
            f.close()                


    def PrintLines(self):
        for line in self.code_generation_table:
            line.LineIntrlPrint()
        
    def PrintLines_InitializationList(self):
        for line in self.InitializationList:
            line.LineIntrlPrint()
        
    def readfield(self,tline,bmfId):
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


    def readoperations(self,bmfId):
        #Parses the netlist assignments into a graph (matrix) 

        global code_generation_top,PRINT_WARNING 
        
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
                outIdx = self.varLegendCol.index(out)
                input_idxs = []
                for k in range(0,num_of_op - 1):
                    inIdx = self.varLegendRow.index(operands[k])
                    self.NodesList[inIdx].AddOutEdge(outIdx + self.i,1) #TODO - remove self.i
                    self.NodesList[outIdx + self.i].AddInEdge(inIdx,1)
                    input_idxs.append(inIdx) #Gate inputs list
                self.Insert_readoperations_parameters(outIdx + (self.lr -self.lc),input_idxs,op) #For statistics      
            tline = bmfId.readline()            

    def GetRoots_list(self): 
        #Returns the graph roots. Also calculates the FO array.
        roots = []
        for i,node in enumerate(self.NodesList):
            row_sum = node.GetNumOfOutEdges()
            if (row_sum == 0):
                if ((node.GetNumOfInEdges()) != 0):
                    roots.append(i)
                else:
                    if (PRINT_WARNING == True):
                        print('** Warning **',self.varLegendRow[i],'has no input')
                    self.Increase_NoInputWireNum_by_one()
                    self.Inset_To_NoInputWireList(i)
            else:
                self.NodesList[i].SetNodeFO(row_sum) 
        print('\n')
        return roots    

    def GetParents_list(self,V_i):
        return self.NodesList[V_i].GetOutEdgesList(True)

    def GetChildrens_list(self,V_i): 
        return self.NodesList[V_i].GetInEdgesList(True)
    
    def ChildrenWithoutInputs_list(self,V_i):
        #Returns the childrens without netlist inputs
    
        childrens = self.GetChildrens_list(V_i)
        childrens_without_inputs = [child for child in childrens if (child in self.LEAFS_inputs) == False]
        return childrens_without_inputs

    def computeCU(self,V_i):
        #Computes the cell usage (CU) of gate (node) V_i  
        if (self.NodesList[V_i].GetNodeCu() > 0):
            return # CU[V_i] was already generated and therefore doesn't change   
        childrens = self.ChildrenWithoutInputs_list(V_i)
        if(len(childrens) == 0): #V_i has no childrens -> V_i is connected to function inputs only
            self.NodesList[V_i].SetNodeCu(1)
        else:
            if (len(childrens) == 1):
                self.computeCU(childrens[0])
                tmp = self.NodesList[childrens[0]].GetNodeCu()
                self.NodesList[V_i].SetNodeCu(tmp)
            else:
                childrens_cu = []
                for child in childrens:
                    self.computeCU(child)
                    childrens_cu.append(self.NodesList[child].GetNodeCu())
                childrens_cu.sort(key=None, reverse=True)               
                num_of_childrens = range(0,len(childrens)) # equal to + (i - 1) for all i in 1 to N (N is the number of childrens)           
                self.NodesList[V_i].SetNodeCu(max(np.add(childrens_cu,num_of_childrens)))
            

    def AllocateRow(self,V_i):
        #Allocates cells to the gate V_i and his children (a sub-tree rooted by V_i). 
        #In a case the allocation for one of V_i's children or V_i itself is failed, the function returns False. On successful allocation returns True.
        childrens = self.ChildrenWithoutInputs_list(V_i) #Equal to C(V_i) - the set of V_i's childrens
        childrens_sorted_by_cu = [[child,self.NodesList[child].GetNodeCu()] for child in childrens] # the loop creates a list composed of pairs of the form [child number, CU[child number]]
        childrens_sorted_by_cu.sort(key = lambda k: k[1], reverse=True) #sorting by CU   
        childrens_sorted_by_cu = [elm[0] for elm in childrens_sorted_by_cu] #taking only the child (vertex) number       
        for V_j in childrens_sorted_by_cu:            
            if (self.NodesList[V_j].GetNodeMap() == 0):
                if (self.AllocateRow(V_j) == False):
                    return False
        if (self.NodesList[V_i].GetNodeMap() == 0): #V_i is not mapped
            self.NodesList[V_i].SetNodeMap(self.AllocateCell(V_i))    
            if (self.NodesList[V_i].GetNodeMap() == 0): #V_i could not be mapped
                return False
        return True

    def AllocateCell(self,V_i):
        FreeCell = self.cells.GetFirst_Available()
        if FreeCell == None:
            if self.cells.IsNotEmpty_Init():
                self.t += 1
                self.Add_To_Initialization_List(self.t, self.cells.init_list_for_json)
                self.cells.Concatenate_init_to_available_list()
                self.cells.Empty_Init()
                FreeCell = self.cells.GetFirst_Available()
            else:
                return 0
        self.cells.DeleteFirst_Available()
        self.cells.Insert_Used(FreeCell,V_i)
        self.t += 1
        self.Insert_AllocateCell_parameters(V_i,FreeCell)
        for V_k in self.ChildrenWithoutInputs_list(V_i):
            self.NodesList[V_k].SetNodeFO(self.NodesList[V_k].GetNodeFO() - 1)
            if (self.NodesList[V_k].GetNodeFO() == 0):
                cell_to_be_moved = self.NodesList[V_k].GetNodeMap()
                self.cells.Delete_Used(cell_to_be_moved)
                self.cells.Insert_Init(cell_to_be_moved)
        return FreeCell

              
    def IncreaseOutputsFo(self):
        # This function created to make sure outputs cells will not evacuated.
        #It is done by increasing their FO by 1

        for idx in range(self.len_input_and_wire,self.lr): #outputs idx range
            self.NodesList[idx].SetNodeFO(self.NodesList[idx].GetNodeFO() + 1)
            
            
    def RunAlgorithm(self):
            #================ SIMPLER algorithm Starts ================ 
            
            #FO array initialized in __init__
            #Cell Usage array initialized in __init__
            #Map is the number of the cell/column V_i is mapped to. Array initialized in __init__      
            #CELLS list initialization 
            ROOTs = self.GetRoots_list() #set of all roots of the graph. Also calculates FO values.
            self.IncreaseOutputsFo() # To ensure outputs who are also inputs, will not be evacuated
            self.t = 0 #Number of clock cycles        

            for netlist_input in range(0,self.i):
                self.cells.Insert_Used(netlist_input,netlist_input)
            for cell in range(self.i,self.N):
                self.cells.Insert_Available(cell) 
            #alg start here
            t1 = time.time()#time
            for r in ROOTs:
                self.computeCU(r)
            if SORT_ROOTS == 'NO':
                for r in ROOTs:    
                    if (self.AllocateRow(r) == False):
                        print('\\\\\\\\\\\\ MAPPING OF',self.Benchmark,'WITH ROW SIZE =',self.N,' \\\\\\\\\\\\\n')
                        print('False - no mapping\n')
                        #code_generation_success_flag = False #Printing flag 
                        #break #To enable multiple runs. To fit the code to the article, comment this line, and uncomment the two next lines. 
                        return False #Cannot find mapping
                t2 = time.time()
                print('time is:',t2-t1)#time
                return True #A mapping of the entire netlist was found 
            else:        
                sorted_ROOTs = [[r,self.NodesList[r].GetNodeCu()] for r in ROOTs]
                if SORT_ROOTS == 'DESCEND':
                    sorted_ROOTs.sort(key = lambda k: k[1], reverse=True)
                elif SORT_ROOTS == 'ASCEND':
                    sorted_ROOTs.sort(key = lambda k: k[1], reverse=False)
                for sr in sorted_ROOTs:    
                    if (self.AllocateRow(sr[0]) == False):
                        print('\\\\\\\\\\\\ MAPPING OF',self.Benchmark,'WITH ROW SIZE =',self.N,' \\\\\\\\\\\\\n')
                        print('False - no mapping\n')
                        return False #cannot find mapping
                return True          


#============ End of Globals variables and Classes ==============


   
#======================== SIMPLER MAPPING =======================
def SIMPLER_Main (BenchmarkStrings, Max_num_gates, ROW_SIZE, Benchmark_name, generate_json, print_mapping, print_warnings):
    global JSON_CODE_GEN, PRINT_CODE_GEN, PRINT_WARNING, SORT_ROOTS
    
    #print controls
    JSON_CODE_GEN = generate_json
    PRINT_CODE_GEN = print_mapping
    PRINT_WARNING = print_warnings
    SORT_ROOTS = 'NO' #Set to one of the follows: 'NO, 'ASCEND' 'DESCEND' 
    
    for Row_size in ROW_SIZE: 
        for Benchmark in BenchmarkStrings:
            
            #Parse operations 
            bmfId = open(Benchmark,"r") #open file       
            SIMPLER_TDS = SIMPLER_Top_Data_Structure(Row_size,bmfId,Benchmark_name)
                          
            if (SIMPLER_TDS.Get_lr()>Max_num_gates or SIMPLER_TDS.Get_lc()>Max_num_gates):
                print("** net too big, skip " + str(SIMPLER_TDS.Get_lr()) +" X " + str(SIMPLER_TDS.Get_lc()) + "\n")
                continue
                              
            #Statistics calculations 
            SIMPLER_TDS.Set_Max_Num_Of_Used_Cells(CellInfo.get_max_num_of_used_cells())
            code_generation_success_flag =SIMPLER_TDS.RunAlgorithm()
            if (code_generation_success_flag == True):
                SIMPLER_TDS.PrintCodeGeneration() 
            
            #Benchmark's end 
            bmfId.close() #close file
            CellInfo.Set_cur_num_of_used_cells_to_zero() #need to initiate because its a class variable
            CellInfo.Set_max_num_of_used_cells_to_zero() #need to initiate because its a class variable
            print('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ \n')

#=========================== End of SIMPLER MAPPING ===========================

#============================== End of code ===================================
