import numpy as np
filename = './BOITEDON/JPBCFIG.DON'

with open(filename, "r") as ins:
    array = []
    for line in ins:
        array.append(line)


#print(array)

import re
re_dbl_fort = re.compile(r'(\d*\.\d+)[dD]([-+]?\d+)')

text = 'DEW=[0.242D+03 -4.320D-06]'
re_dbl_fort.sub(r'\1e\2', text)


def get_next_integer(line_list, current):

    for i in range(current,len(line_list)):
        try :
            print(line_list[i])
            a = int(line_list[i])
            success =True
        except:
            success =False
        if success:
            return a , i

def get_next_double(line_list, current):

    for i in range(current,len(line_list)):
        try :
            print(line_list[i])
            a = float(re_dbl_fort.sub(r'\1e\2', line_list[i]))
            success =True
        except:
            success =False
        if success:
            return a , i


elem = []

l=0
for line in array:
    if (line[0:5] == 'ELMxx'):
        number_elem, current = get_next_integer(array[l+1].split(' '), 0) 
        print('number_elem', number_elem, current)
        input()
        elem_start= l+3
        elem_end= l+3+number_elem
        print('number_elem', number_elem)
    if (line[0:5] == 'IMMxx'):
        print(line)
        number_node, current = get_next_integer(array[l+1].split(' '), 0) 
        node_start= l+3
        node_end= l+3+number_node
        print('node_start',node_start)
        print('node_end',node_end)
        print('number_node', number_node)      
    l =l+1

input()
for i in range(elem_start,elem_end,1):
    elem_info=array[i].split(' ')
    print(elem_info)
    elem_current= {}
    elem_nodes = []
    elem_type=elem_info[-5][0:5]
    print('elem_type',elem_type)
    elem_mat=elem_info[-3][0:5]
    print('elem_mat',elem_mat)
    elem_bdy=elem_info[-1][0:5]
    print('elem_bdy',elem_bdy)

    if (elem_type == 'T3xxx'):
        for a in elem_info:
            try :
                elem_nodes.append(int(a)-1)
            except:
                pass
        elem_current['nodes'] = elem_nodes[1:]
        elem_current['type'] = elem_type
        elem_current['mat'] = elem_mat
        elem_current['bdy'] = elem_bdy
        elem_current['number'] = elem_nodes[0]-1
        
        elem.append(elem_current)

node = []



for i in range(node_start,node_end,1):
    node_info=array[i].split(' ')
    print('node_info', node_info)
    node_current = {}
    node_coord = []
    
    node_type=node_info[-1][0:5]
    print('node_type',node_type)
    n=0

    if (node_type == 'NO2xx'):
        x , current =get_next_double(node_info, 0)
        y , current =get_next_double(node_info, current+1)

        node_coord.append(x)
        node_coord.append(y)
        ii , current =get_next_integer(node_info, current+1)
        node_coord.append(ii-1)
        print(x,y, ii-1)
        
        print(node_coord)
        node_current['number']= node_coord[-1]
        node_current['type'] =node_type
        node_current['coord']=node_coord[0:2]
        node.append(node_current)
        
    #input()


  


print(elem)
print(node)
import pickle
pickle.dump(elem, open('elem.dat', 'wb'))
pickle.dump(node, open('node.dat', 'wb'))
