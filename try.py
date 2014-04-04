from sets import Set

import networkx as nx;
'''
for line in open("users.txt", 'r'):
	print len(line.split(','));

for line in open("users.txt", 'r'):
	print len(line.split(','));

'''


tot_users = 0;
G = nx.Graph();
users = dict();
for line in open("users.txt", 'r'):
	for i in line.split(','):
		users[i] = tot_users;
		G.add_node(tot_users);
		tot_users += 1;
print tot_users;
for line in open("edges.csv",'r'):
	tem = line.split();
	if (tem[0] in users) and (tem[1] in users):
		G.add_edge(users[tem[0]], users[tem[1]]);

print "Number of Connected Components:", nx.number_connected_components(G);
cc = nx.connected_components(G)
print nx.graph_clique_number(G);

exit();

