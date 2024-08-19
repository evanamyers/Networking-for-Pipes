import sys
import networkx as nx
import matplotlib.pyplot as plt
import pickle

# # Save the graph to a file
with open('graph.pkl', 'rb') as f:
    G = pickle.load(f)

fig, ax = plt.subplots()

# Draw the edges
for edge in G.edges():
    ax.plot(*zip(edge[0], edge[1]), color='gray')

# Draw the nodes
for node in G.nodes():
    ax.plot(*node, 'o', color='lightblue')
    #ax.text(node[0], node[1], f'({node[0]:.2f}, {node[1]:.2f})', ha='right', va='top', fontsize=8)

# Show the plot
plt.show()



# fig, ax = plt.subplots()
# pos = {node_id: coords for node_id, coords in G.nodes(data='coords')}
# nx.draw(loaded_graph, pos=pos, with_labels=False, ax=ax)
# plt.show()

# pos = nx.spring_layout(loaded_graph)
# pos = nx.get_node_attributes(loaded_graph, 'pos')  # Get the node coordinates

# nx.draw(G, pos=pos, arrows=True, node_size=500, node_color='lightblue', edge_color='gray')
# plt.show()

# nx.spring_layout(loaded_graph)
# plt.show()

