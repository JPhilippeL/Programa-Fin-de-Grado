class GraphManager:
    @staticmethod
    def delete_node(graph, node_id):
        if node_id in graph:
            graph.remove_node(node_id)

    @staticmethod
    def modify_node(graph, node_id, new_element):
        if node_id in graph:
            graph.nodes[node_id]['element'] = new_element
