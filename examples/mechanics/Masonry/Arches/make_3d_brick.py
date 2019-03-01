import pickle
import numpy as np
elem = list(pickle.load(open( 'elem.dat', "rb" )))
node = list(pickle.load(open( 'node.dat', "rb" )))
print(elem[0])
# for e in elem:
#     print(e)

# compute coordination
n_node=0
for n in node:
    n_node =n_node+1
    
nodes_vector = np.zeros((n_node+1,1), dtype='int')

nodes_elements = []

for i in range(n_node+1):
    nodes_elements.append( [])
#print(nodes_vector)

#print(elem)


for e in elem:
    print(e)
    nodes_e = e['nodes']
    print('nodes_e', nodes_e)
    
    for n in nodes_e:
        nodes_vector[n] =nodes_vector[n]+1
        nodes_elements[n].append(e)
    #input()
         
#print(nodes_vector)
#print(nodes_elements)

brick = []
for n_e in range(n_node):
    if nodes_vector[n_e] ==8 :
        #print('n_e',n_e)
        element_nodes = nodes_elements[n_e]
        vertices = []
        for e in element_nodes:
            print('e[nodes]', e['nodes'])
            vertices.extend(e['nodes'])
            if n_e in vertices:
                vertices.remove(n_e)
        vertices= list(set(vertices))
        brick.append(vertices)
        #print(n_e)

print(brick)

n_n = n_node

new_nodes = []


new_brick = []

i=0 
for b in brick:
    v=[]
    print(b)
    for n in b:
        new_node0={}
        new_node0['number']=n
        #input()
        coord = list(node[n]['coord']) #weird
        new_node0['coord'] = [coord[0], 0.0, coord[1]  ]
        v.append(n)
        new_nodes.append(new_node0)
        print('new_node0',new_node0)
    for n in b:
        n_n =n_n+1      
        new_node1={}
        new_node1['number']=n_n
        coord = list(node[n]['coord']) #weird
        new_node1['coord'] = [ coord[0], 0.25, coord[1]  ]
        new_nodes.append(new_node1)
        print('new_node1',new_node1)
        
        v.append(n_n)
    print('v',v)
    new_brick.append(v)
    i=i+1
    # if (i > 10 ):
    #     break
    

    #input()

print(new_nodes)       
print(new_brick)       

pickle.dump(new_nodes, open( 'vertices.dat', "wb" ))
pickle.dump(new_brick, open( 'brick.dat', "wb" ))

