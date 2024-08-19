

def createGraph():

    diameterField = 'DIAMETER'
    materialField = 'MATERIAL'
    fields = ["GlobalID", "SHAPE@", diameterField, materialField]
    lineNetworkDict = {}
    lineCoordList = []
    lineShapeList = []
    for l, each in enumerate(listLineFC):
        letter = chr(ord('a') + l)
        i = 0
        with arcpy.da.SearchCursor(each, fields, where_clause=lineFilter) as cur:
            for row in cur:
                lineXY = []
                nodes = []
                for part in row[1]:
                    for pnt in part:
                        coord = (float("{:.10f}".format(pnt.X)), float("{:.10f}".format(pnt.Y)))
                        lineXY.append(coord)
                        nodeID = f'{letter}{i}'
                        # G.add_node(nodeID, coords=coord)
                        nodes.append(nodeID)
                        lineCoordList.append(coord)
                        i += 1
                lineShapeList.append(lineXY)
                lineNetworkDict[row[0]] = (nodes, lineXY, row[2], row[3])

    arcpy.AddMessage('Collected line features')
    # Add nodes and edges
    for node_id, (labels, coordinates, diameter, material) in lineNetworkDict.items():
        for label, coord in zip(labels, coordinates):
            G.add_node(coord, id=node_id, label=label)

        for i in range(len(coordinates) - 1):
            x1, y1 = coordinates[i]
            x2, y2 = coordinates[i + 1]
            # add weights to each edge based on length
            G.add_edge(coordinates[i], coordinates[i + 1], weight=math.sqrt((x2 - x1)**2 + (y2 - y1)**2),
                       diameter=diameter, material=material)

    arcpy.AddMessage('Created nodes and edges')
    if G:

        edgesCoords = [edge for edge in G.edges()]

        if edgesCoords:
            output_fcEdges = arcpy.CreateFeatureclass_management(arcpy.env.workspace, 'network_edges', 'POLYLINE',
                                                                 spatial_reference=arcpy.Describe(
                                                                     listLineFC[0]).SpatialReference)
            arcpy.management.AddFields(in_table='network_edges',
                                       field_description="Diameter DOUBLE Diameter 10 # # ;Material TEXT Material 255 # #")
            dirname = os.path.dirname(arcpy.Describe(output_fcEdges).catalogPath)
            networkFields = ["SHAPE@", diameterField, materialField]
            with arcpy.da.Editor(dirname):
                with arcpy.da.InsertCursor(output_fcEdges, networkFields) as cursor:  # ['SHAPE@', diameterField]
                    for coords, values in G.edges.items():
                        line = LineString(coords)
                        polyline = arcpy.FromWKT(line.wkt)  # Unpack the tuple
                        cursor.insertRow([polyline, values.get(diameterField.lower()), values.get(materialField.lower())])

            currentMap.addDataFromPath(output_fcEdges)
        arcpy.AddMessage('Created Graph')

    return G


def fixBrokenConnections():
    # try:
    print("Graph is not connected.")

    # Get the connected components
    if multiDirection == 1:
        # connected_components = list(nx.connected_components(G))
        largest_component = max(nx.connected_components(G), key=len)  # this lists just the nodes
        largest_component_edges = G.subgraph(largest_component).edges()  # this lists just the edges

        # Find nodes that are not part of the largest connected component
        isolated_nodesCoords = [node for node in G.nodes() if node not in largest_component]
        isolated_edgesCoords = [edge for edge in G.edges() if edge not in largest_component_edges]
    if multiDirection == 0:
        largest_scc = max(nx.strongly_connected_components(G), key=len)
        largest_wcc = max(nx.weakly_connected_components(G), key=len)
        largest_component_edges = G.subgraph(largest_wcc).edges()

        isolated_nodesCoords = [node for node in G.nodes() if node not in largest_wcc]
        isolated_edgesCoords = [edge for edge in G.edges() if edge not in largest_component_edges]


    # Add isolated nodes to the map if there are any
    # if isolated_nodesCoords:
    #     output_fcNodes = arcpy.CreateFeatureclass_management(arcpy.env.workspace, 'isolated_nodes', 'POINT',
    #                                                          spatial_reference=arcpy.Describe(
    #                                                              listLineFC[0]).SpatialReference)
    #     dirname = os.path.dirname(arcpy.Describe(output_fcNodes).catalogPath)
    #     with arcpy.da.Editor(dirname, multiuser_mode=True):
    #         with arcpy.da.InsertCursor(output_fcNodes, ['SHAPE@']) as cursor:
    #             for coords in isolated_nodesCoords:
    #                 point = arcpy.Point(*coords)  # Unpack the tuple
    #                 cursor.insertRow([point])
    #
    #     currentMap.addDataFromPath(output_fcNodes)

    # Add isolated edges to the map if there are any
    if isolated_edgesCoords:
        output_fcEdges = arcpy.CreateFeatureclass_management(arcpy.env.workspace, 'isolated_lines', 'POLYLINE',
                                                             spatial_reference=arcpy.Describe(
                                                                 listLineFC[0]).SpatialReference)
        dirname = os.path.dirname(arcpy.Describe(output_fcEdges).catalogPath)
        with arcpy.da.Editor(dirname, multiuser_mode=True):
            with arcpy.da.InsertCursor(output_fcEdges, ['SHAPE@']) as cursor:
                for coords in isolated_edgesCoords:
                    line = LineString(coords)
                    polyline = arcpy.FromWKT(line.wkt)  # create an arcgis geometry
                    cursor.insertRow([polyline])

        currentMap.addDataFromPath(output_fcEdges)




if __name__ == '__main__':

    import arcpy
    import os
    import math
    import networkx as nx
    import pickle
    from shapely.geometry import Point, LineString

    multiDirection = arcpy.GetParameter(0)
    aprx = arcpy.mp.ArcGISProject('current')
    currentMap = aprx.activeMap
    arcpy.env.workspace = aprx.defaultGeodatabase

    lineNetworkFeatures = arcpy.GetParameterAsText(1).replace("'", "")
    listLineFC = lineNetworkFeatures.split(";") if ";" in lineNetworkFeatures else [lineNetworkFeatures]
    lineFilter = arcpy.GetParameter(2)  # Definition queries do not use "GetParameterAsText"

    arcpy.AddMessage(lineFilter)
    arcpy.AddMessage(str(listLineFC))

    junctions = arcpy.GetParameterAsText(3).replace("'", "")  # Should be a Fitting feature
    junctionsFC = junctions.split(";") if ";" in junctions else [junctions]
    junctionFilter = arcpy.GetParameter(4)


    # Create a graph of line's nodes (vertex)
    if multiDirection == 1:
        arcpy.AddMessage("Creating Multi-Directional Network")
        print("Creating Multi-Directional Network")
        G = nx.Graph()  # used to be MultiGraph
        G.graph["crs"] = "EPSG:2236"
        createGraph()

        isConnected = nx.is_connected(G)
        arcpy.AddMessage(f"Is connected: {isConnected}")
        print("Is connected:", isConnected)
        if not isConnected:
            arcpy.AddMessage('Some parts of Network are disconnected.  Creating an "Isolated Lines" layer to show areas not '
                             'connected to larger network.')
            fixBrokenConnections()

    if multiDirection == 0:
        arcpy.AddMessage("Creating Single Direction Network")
        print("Creating Single Direction Network")
        G = nx.DiGraph()  # Di meaning 'Directed' (only one direction)  # used to be MultiDiGraph
        G.graph["crs"] = "EPSG:2236"
        createGraph()

        # NOTE:  strongly_connected is not used because there are naturally dead ends in a pipe network
        # is_strongly_connected = nx.is_strongly_connected(G)
        is_weakly_connected = nx.is_weakly_connected(G)
        # arcpy.AddMessage(f"Strongly connected: {is_strongly_connected}")
        arcpy.AddMessage(f"Weakly connected: {is_weakly_connected}")
        # print(f"Strongly connected: {is_strongly_connected}")
        print(f"Weakly connected: {is_weakly_connected}")
        if not is_weakly_connected:
            arcpy.AddMessage('Some parts of Network are disconnected.  Creating an "Isolated Lines" layer to show areas not '
                             'connected to larger network.')
            fixBrokenConnections()

    print(G.edges())
    print(G.nodes())

    # Save the graph to a file
    with open('graph.pkl', 'wb') as f:
        pickle.dump(G, f)



