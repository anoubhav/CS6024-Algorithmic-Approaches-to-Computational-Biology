import sys
import os
# os.system('cls')

def main():

    # O(n) time complexity, where n is the number of input k-mers
    edges = list()
    for k_mer in sys.stdin:
        if k_mer == '\n': break
        k_mer = k_mer.strip()
        edges.append((k_mer[:-1], k_mer[1:]))
    
    # O((n log n)*2*(k-1)) i.e. O(k*nlogn) time complexity for sorting into lexicographical order
    edges = sorted(edges, key=lambda x: x[0]+x[1])
    
    # O(n) time for printing
    curr_node = edges[0][0] # stores the current node
    adj_nodes = list()
    tot_len   = len(edges)
    index     = 0

    for node, adj_node in edges:
        if node == curr_node and index != tot_len-1:
            adj_nodes.append(adj_node)
        elif index == tot_len - 1:
            if node == curr_node:
                adj_nodes.append(adj_node)
                print('{} -> {}'.format(curr_node, ','.join(adj_nodes)))
            else:
                print('{} -> {}'.format(curr_node, ','.join(adj_nodes)))
                print('{} -> {}'.format(node, adj_node))

        else:
            print('{} -> {}'.format(curr_node, ','.join(adj_nodes)))
            curr_node = node
            adj_nodes = list()
            adj_nodes.append(adj_node)
        index += 1

if __name__ == "__main__":
    main()