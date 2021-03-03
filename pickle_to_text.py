import pickle

c = pickle.load(open('','rb'))
cd = c.to_node_community_map()

with open('out.txt','w') as outfile:
    for key in cd.keys():
        outfile.write(key + '\t' + '\t'.join(cd[key]))

