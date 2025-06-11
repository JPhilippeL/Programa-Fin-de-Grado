import networkx as nx
class GraphManager:
    @staticmethod
    def delete_node(graph, node_id):
        if node_id in graph:
            graph.remove_node(node_id)

    @staticmethod
    def modify_node(graph, node_id, new_element):
        if node_id in graph:
            graph.nodes[node_id]['element'] = new_element

    @staticmethod
    def add_edge(graph, source_id, target_id, bond_type="SINGLE"):
        if not graph.has_edge(source_id, target_id):
            graph.add_edge(source_id, target_id, bond_type=bond_type)
