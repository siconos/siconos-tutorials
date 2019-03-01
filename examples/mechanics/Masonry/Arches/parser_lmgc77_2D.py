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



elem = []

l=0
for line in array:
    if (line[0:5] == 'ELMxx'):
        print(line)
        elem_info=array[l+1].split(' ')

        number_elem = int(elem_info[2])
        elem_start= l+3
        elem_end= l+3+number_elem
        print('number_elem', number_elem)
    if (line[0:5] == 'IMMxx'):
        print(line)
        node_info=array[l+1].split(' ')
        number_node = int(node_info[2])
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
    if (node_type == 'NO2xx'):
        node_coord.append(int(node_info[-3]))
    n=0
    for a in node_info:
        #print(a)
        #print(re_dbl_fort.sub(r'\1e\2', a))
        #input()
        if (node_type == 'NO2xx'):
            try :
                node_coord.append(float(re_dbl_fort.sub(r'\1e\2', a)))
                succeed=True
            except:
                succeed = False
                pass
            if succeed :
                n=n+1
        if (n == 2):
            #print('good number of nodes')
            break

            
   
    #input()

    if (node_type == 'NO2xx'):
        print(node_coord)
        node_current['number']= node_coord[0]-1
        node_current['type'] =node_type
        node_current['coord']=node_coord[1:3]
        node.append(node_current)
        #input()


print(elem)
print(node)
import pickle
pickle.dump(elem, open('elem.dat', 'wb'))
pickle.dump(node, open('node.dat', 'wb'))
