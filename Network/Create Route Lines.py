import networkx as nx
from networkx.exception import NetworkXNoPath, NodeNotFound
# import arcpy
# import os
# import networkx as nx
# import pickle
# from shapely.geometry import Point, LineString
# from shapely import STRtree
# from shapely.ops import split
#
#
# arcpy.SetLogMetadata(False)
# arcpy.SetLogHistory(False)
# arcpy.env.overwriteOutput = True
#
startFeature = "wControlValve"
startIDField = 'MXLOCATION'

endFeature = "wSystemValve"
endIDField = 'MXLOCATION'

startPointFilter = None
endPointFilter = None

aprx = arcpy.mp.ArcGISProject('current')
currentMap = aprx.activeMap
arcpy.env.workspace = aprx.defaultGeodatabase
ws = arcpy.env.workspace
graphFile = r"W:\GIS\WORKING FILES\EM Working Files\All Files for Projects\Python\Network\graph.pkl"
with open(graphFile, 'rb') as f:
    G = pickle.load(f)


def getLocations():

    # Create a tree for the Graph Nodes
    graphNodeList = [Point(node[0]) for node in G.nodes(data=True)]
    graphNodeTree = STRtree(graphNodeList)

    # Get starting point coordinate:
    startPointDict = {}
    startPointList = []
    # soidField = arcpy.Describe(startFeature).OIDFieldName
    with arcpy.da.SearchCursor(startFeature, [startIDField, "SHAPE@XY"], where_clause=startPointFilter) as cursor:
        for row in cursor:
            if row[0]:
                x, y = row[1]
                coords = (float("{:.10f}".format(x)), float("{:.10f}".format(y)))
                startPointList.append(Point(coords))
                startPointDict[row[0]] = coords

    startPointTree = STRtree(startPointList)

    # Get end point coordinate:
    endPointDict = {}
    endPointList = []
    # eoidField = arcpy.Describe(endFeature).OIDFieldName
    with arcpy.da.SearchCursor(endFeature, [endIDField, "SHAPE@XY"], where_clause=endPointFilter) as cursor:
        for row in cursor:
            if row[0]:
                x, y = row[1]
                coords = (float("{:.10f}".format(x)), float("{:.10f}".format(y)))
                endPointList.append(Point(coords))
                endPointDict[row[0]] = coords

    endPointTree = STRtree(endPointList)

    try:
        startNodeList = []  # Shapely Point objects
        startNodeDict = {}  # MXLOCATION and coordinate of node
        for key, values in startPointDict.items():
            indices = graphNodeTree.query(Point(values))
            if indices:
                intersecting_nodes = [graphNodeList[i] for i in indices][0]
                if intersecting_nodes:
                    intersectPoint = Point(intersecting_nodes.coords[0])
                    startNodeList.append(intersectPoint)
                    startNodeDict[key] = intersecting_nodes.coords[0]
        startNodeTree = STRtree(startNodeList)
    except IndexError:
        raise "The Start Feature was not found on a vertex. Please snap the feature to a vertex on a related line."
    try:
        endNodeList = []  # Shapely Point objects
        endNodeDict = {}  # MXLOCATION and coordinate of node
        for key, values in endPointDict.items():
            indices = graphNodeTree.query(Point(values))
            if indices:
                intersecting_nodes = [graphNodeList[i] for i in indices][0]
                if intersecting_nodes:
                    intersectPoint = Point(intersecting_nodes.coords[0])
                    endNodeList.append(intersectPoint)
                    # for each in startNodeList:
                    #     dist = each.distance(intersectPoint)
                    endNodeDict[key] = intersecting_nodes.coords[0], 0  # , dist

        # endNodeDict = {k: v for k, v in sorted(endNodeDict.items(), key=lambda item: item[1][1])}
        # endNodeTree = STRtree(endNodeList)
    except IndexError:
        raise "The End Feature was not found on a vertex. Please snap the feature to a vertex on a related line."

    # print(startNodeList)

    return startNodeDict, endNodeDict, graphNodeList, graphNodeTree


def createShortestPath(startNodeDict, endNodeDict, graphNodeList, graphNodeTree):

    # Create list to store the coordinates of paths
    paths = []
    pathNum = 1
    pathsDict = {}
    for startkey, startvalues in startNodeDict.items():
        for endkey, endvalues in endNodeDict.items():
            try:
                indices = graphNodeTree.query(Point(endvalues[0]))
                if indices:
                    intersecting_nodes = [graphNodeList[i] for i in indices][0]
                    if intersecting_nodes:
                        intersectPoint = Point(intersecting_nodes.coords[0])
                        dist = Point(startvalues).distance(intersectPoint)
                        endNodeDict[endkey] = intersecting_nodes.coords[0], dist
            # endNodeDict = {k: v for k, v in sorted(endNodeDict.items(), key=lambda item: item[1][1])}
                shortestPath = nx.shortest_path(G, startvalues, endvalues[0], weight='weight')  # , weight='weight'
                pathCoords = [(float("{:.10f}".format(x)), float("{:.10f}".format(y))) for x, y in shortestPath]
                pathShape = LineString(pathCoords)
                paths.append(pathShape)  # extend
                pathsDict[pathNum] = pathShape, pathShape.length, startkey, endkey
                pathNum += 1
            except (NetworkXNoPath, NodeNotFound):
                print(f"No route from {startkey} to {endkey}. Skipping")
                pass

    # try:
    #     for startkey, startvalues in startNodeDict.items():
    #         for endkey, endvalues in endNodeDict.items():
    #             allPaths = nx.all_simple_paths(G, startvalues, endvalues)
    #     # return shortestPath
    # except():
    #     print("This network is still broken.")

    return pathsDict


def isolationPath(startNodeDict, endNodeDict, graphNodeList, graphNodeTree):

    arcpy.AddMessage("Finished searching paths")

    pathNum = 1
    pathsDict = {}
    for startkey, startvalue in startNodeDict.items():
        attempts = 1
        print(attempts)
        finaledgeList = []
        for endkey, endvalue in endNodeDict.items():
            try:
                path = nx.shortest_path(G, startvalue, endvalue[0], weight='weight')  # (G, source, target)
                pathShape = LineString(path)  # remove if linestring doesn't work
                if any(point in list(values[0] for values in endNodeDict.values()) for point in path[1:-1]):
                    pass
                else:
                    pathsDict[pathNum] = pathShape, pathShape.length, startkey, endkey
                    pathNum += 1
                    # Get network edges that touch the path
                    edgeList = [(u, v) for u, v in G.edges(path)]
                    # Remove network edges that are underneath the path
                    edges_to_remove = [edge for edge in edgeList if any(
                        edge[0] in pathShape.coords and edge[1] in pathShape.coords for pathShape, _, _, _ in
                        pathsDict.values())]
                    for edge in edges_to_remove:
                        edgeList.remove(edge)
                    availRoutes = []
                    # Remove network edges that touch the endpoint of the path, leaving only open available routes
                    for u, v in edgeList[:]:
                        if endvalue[0] in [u, v]:
                            edgeList.remove((u, v))
                        else:
                            edge_line = LineString([u, v])
                            availRoutes.append(edge_line)
                    for each in edgeList:
                        if each not in finaledgeList:
                            finaledgeList.append(each)
                # Keep track of the edges that haven't been used in a path, using 'finaledgeList'
                for value in pathsDict.values():
                    line = value[0]  # Get the LineString object
                    valueCoords = [(x, y) for x, y in line.coords]  # Extract the coordinates as a list of tuples
                    found = False
                    for edge in finaledgeList:
                        if any((valueCoords[i], valueCoords[i + 1]) == edge for i in range(len(valueCoords) - 1)):
                            found = True
                            break
                    if found:
                        finaledgeList.remove(edge)
                    else:
                        continue
                    # if len(finaledgeList) == 0:
                    #     print('fully isolated!@!@@')
                    #     break
                    # else:
                    #     continue
                # if finaledge list is at 1
                if len(finaledgeList) >= 1:
                    attempts += 1
                    print(attempts)
                if len(finaledgeList) == 0:
                    print('fully isolated!@!@!@!')
                    break
                if attempts == 1000:
                    print('could not fully isolate. check route layer for open ends.')
                    break
            # except nx.exception.NetworkXNoPath or networkx.exception.NodeNotFound:
            except (NetworkXNoPath, NodeNotFound) as e:
                print(f"Skipping due to error: {e}")
                pass

    return pathsDict


def createPathLayer(pathsDict):

    # Create a new polyline feature class
    arcpy.AddMessage("Creating Route layer")
    # arcpy.env.addOutputsToMap = True
    pathFC = arcpy.CreateFeatureclass_management(ws, 'Route', 'POLYLINE', spatial_reference=arcpy.Describe(startFeature).SpatialReference)

    # pathLyr = arcpy.mp.ArcGISProject('current').activeMap.listLayers('Route')[0]
    # sym = pathLyr.symbology
    # sym.updateRenderer('SimpleRenderer')

    allPathdesc = arcpy.Describe(pathFC)
    dirname = os.path.dirname(arcpy.Describe(pathFC).catalogPath)
    with arcpy.da.Editor(dirname, multiuser_mode=allPathdesc.isVersioned):
        with arcpy.da.InsertCursor(pathFC, ['SHAPE@']) as cursor:
            for item in pathsDict.values():  # newPathList
                pathLine = item[0]
                # turn shapely LineString into an arcgis geometry
                arcgis_geom = arcpy.FromWKT(pathLine.wkt)
                # then add it to feature class
                cursor.insertRow([arcgis_geom])


    # if 'Route' not in currentMap.listLayers():
    currentMap.addDataFromPath(pathFC)


if __name__ == '__main__':

    import arcpy
    import os
    import networkx as nx
    import pickle
    from shapely.geometry import Point, LineString
    from shapely import STRtree
    from shapely.ops import split
    import networkx.exception

    arcpy.SetLogMetadata(False)
    arcpy.SetLogHistory(False)
    arcpy.env.overwriteOutput = True

    isolationRoute = arcpy.GetParameter(0)

    startFeature = arcpy.GetParameterAsText(1)
    startIDField = arcpy.GetParameterAsText(2)

    endFeature = arcpy.GetParameterAsText(3)
    endIDField = arcpy.GetParameterAsText(4)

    startPointFilter = None
    endPointFilter = None
    if startFeature == endFeature:
        startPointFilter = arcpy.GetParameterAsText(5)
        endPointFilter = arcpy.GetParameterAsText(6)

    aprx = arcpy.mp.ArcGISProject('current')
    currentMap = aprx.activeMap
    arcpy.env.workspace = aprx.defaultGeodatabase
    ws = arcpy.env.workspace

    # Load the graph from the file
    with open("graph.pkl", 'rb') as f:
        G = pickle.load(f)

    startNodeDict, endNodeDict, graphNodeList, graphNodeTree = getLocations()

    if startFeature and endFeature:
        if isolationRoute == 1:
            pathsDict = isolationPath(startNodeDict, endNodeDict, graphNodeList, graphNodeTree)
        else:
            pathsDict = createShortestPath(startNodeDict, endNodeDict, graphNodeList, graphNodeTree)

        createPathLayer(pathsDict)

    del aprx, currentMap

