import pickle
import sys
c = pickle.load(open(sys.argv[1],'rb'))
cd = c.to_node_community_map()

with open(sys.argv[2],'w') as outfile:
    for key in cd.keys():
        outfile.write(key + '\t' + '\t'.join([str(x) for x in cd[key]]) + '\n')

