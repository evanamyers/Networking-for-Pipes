

def getLocations():

    # Create a tree for the Graph Nodes
    graphNodeList = [Point(node[0]) for node in G.nodes(data=True)]
    graphNodeTree = STRtree(graphNodeList)

    # Get starting point coordinate:
    startPointDict = {}
    startPointList = []
    soidField = arcpy.Describe(startFeature).OIDFieldName
    with arcpy.da.SearchCursor(startFeature, [soidField, startIDField, "SHAPE@XY"], where_clause=startFeatureFilter) as cursor:
        for row in cursor:
            if row[1]:
                x, y = row[2]
                coords = (float("{:.10f}".format(x)), float("{:.10f}".format(y)))
                startPointList.append(Point(coords))
                startPointDict[row[1]] = coords
            else:
                arcpy.AddWarning(f'A start feature, OID {row[0]}, has an invalid ID value. Skipping')

    # Get end point coordinate:
    endPointDict = {}
    endPointList = []
    eoidField = arcpy.Describe(endFeature).OIDFieldName
    with arcpy.da.SearchCursor(endFeature, [soidField, endIDField, "SHAPE@XY"], where_clause=endFeatureFilter) as cursor:
        for row in cursor:
            if row[1]:
                x, y = row[2]
                coords = (float("{:.10f}".format(x)), float("{:.10f}".format(y)))
                endPointList.append(Point(coords))
                endPointDict[row[1]] = coords
            else:
                arcpy.AddWarning(f'An end feature, OID {row[0]}, has an invalid ID value. Skipping')

    startNodeList = []  # Shapely Point objects
    startNodeDict = {}  # MXLOCATION and coordinate of node
    for key, values in startPointDict.items():
        try:
            indices = graphNodeTree.query(Point(values))
            intersecting_nodes = [graphNodeList[i] for i in indices][0]
            if intersecting_nodes:
                intersectPoint = Point(intersecting_nodes.coords[0])
                startNodeList.append(intersectPoint)
                startNodeDict[key] = intersecting_nodes.coords[0]
        except IndexError:
            print(f"A start feature, {key}, is not connected to the network")
            arcpy.AddWarning(f"A start feature, {key}, is not connected to the network")

    endNodeList = []  # Shapely Point objects
    endNodeDict = {}  # MXLOCATION and coordinate of node
    for key, values in endPointDict.items():
        try:
            indices = graphNodeTree.query(Point(values))
            intersecting_nodes = [graphNodeList[i] for i in indices][0]
            if intersecting_nodes:
                intersectPoint = Point(intersecting_nodes.coords[0])
                endNodeList.append(intersectPoint)
                endNodeDict[key] = intersecting_nodes.coords[0], 0  # , dist
        except IndexError:
            print(f"An end feature, {key}, is not connected to the network")
            arcpy.AddWarning(f"An end feature, {key}, is not connected to the network")

    barrierFeatureDict = {}
    barrierNodeDict = {}
    if barrierFeature:
        barrierFeatureDict = {}
        barrierFeatureList = []
        soidField = arcpy.Describe(startFeature).OIDFieldName
        with arcpy.da.SearchCursor(barrierFeature, [soidField, "SHAPE@XY"], where_clause=barrierFeatureFilter) as cursor:
            for row in cursor:
                if row[0]:
                    x, y = row[1]
                    coords = (float("{:.10f}".format(x)), float("{:.10f}".format(y)))
                    barrierFeatureList.append(Point(coords))
                    barrierFeatureDict[row[0]] = coords

        barrierNodeList = []
        barrierNodeDict = {}
        for key, values in barrierFeatureDict.items():
            try:
                indices = graphNodeTree.query(Point(values))
                intersecting_nodes = [graphNodeList[i] for i in indices][0]
                if intersecting_nodes:
                    intersectPoint = Point(intersecting_nodes.coords[0])
                    barrierNodeList.append(intersectPoint)
                    barrierNodeDict[key] = intersecting_nodes.coords[0]
            except IndexError:
                print(f"A barrier feature, {key}, is not connected to the network")
                arcpy.AddWarning(f"A barrier feature, {key}, is not connected to the network")

    return startNodeDict, endNodeDict, graphNodeList, graphNodeTree, barrierNodeDict


def createShortestPath(startNodeDict, endNodeDict, graphNodeList, graphNodeTree, barrierNodeDict):

    # Create list to store the coordinates of paths
    paths = []
    pathNum = 1
    pathsDict = {}
    for startkey, startvalue in startNodeDict.items():
        # arcpy.AddMessage(f"Searching for {startkey}")
        print(f"Searching for {startkey}")
        # try:
        for endkey, endvalue in endNodeDict.items():
            indices = graphNodeTree.query(Point(endvalue[0]))
            if indices:
                intersecting_nodes = [graphNodeList[i] for i in indices][0]
                if intersecting_nodes:
                    intersectPoint = Point(intersecting_nodes.coords[0])
                    dist = Point(startvalue).distance(intersectPoint)
                    endNodeDict[endkey] = intersecting_nodes.coords[0], dist
        sortedendNodes = sorted(endNodeDict, key=lambda k: endNodeDict[k][1])
        # minKey = min(endNodeDict, key=lambda k: endNodeDict[k][1])
        # print(f'closest end feature {minKey}')
        pathFound = False
        while not pathFound:
            for minKey in sortedendNodes:
                try:
                    # print(endNodeDict[minKey][0])
                    shortestPath = nx.shortest_path(G, startvalue, endNodeDict[minKey][0], weight='weight')
                    # pathCoords = [(float("{:.10f}".format(x)), float("{:.10f}".format(y))) for x, y in shortestPath]
                    pathShape = LineString(shortestPath)
                    if any(point in list(values for values in barrierNodeDict.values()) for point in zip(shortestPath[:-1],
                                                                                               shortestPath[1:])):
                        print("blocked")
                        for node_pair in zip(shortestPath[:-1], shortestPath[1:]):
                            if node_pair in list(values for values in barrierNodeDict.values()):
                                G.remove_edge(*node_pair)
                        continue

                    # arcpy.AddMessage(f"found path from {startkey} to {minKey}")
                    print(f"found path from {startkey} to {minKey}")
                    paths.append(pathShape)  # extend
                    pathsDict[pathNum] = pathShape, pathShape.length, startkey, endkey
                    pathNum += 1
                    pathFound = True
                    break
                except nx.exception.NetworkXNoPath or nx.exception.NodeNotFound:
                    print(f"{startkey} is not part network. Skipping")
                    pass
        if not pathFound:
            print(f"No valid path found for {startkey}.  Check pipe direction or barriers.")
            arcpy.AddWarning(f"No valid path found for {startkey}.  Check pipe direction or barriers.")

    arcpy.AddMessage("Finished searching paths")

    return pathsDict


def isolationPath(startNodeDict, endNodeDict, graphNodeList, graphNodeTree, barrierNodeDict):

    arcpy.AddMessage("Creaing isolation paths.")

    for startkey, startvalue in startNodeDict.items():
        for endkey, endvalue in endNodeDict.items():
            indices = graphNodeTree.query(Point(endvalue[0]))
            if indices:
                intersecting_nodes = [graphNodeList[i] for i in indices][0]
                if intersecting_nodes:
                    intersectPoint = Point(intersecting_nodes.coords[0])
                    dist = Point(startvalue).distance(intersectPoint)
                    endNodeDict[endkey] = intersecting_nodes.coords[0], dist
    endNodeDict = {k: v for k, v in sorted(endNodeDict.items(), key=lambda item: item[1][1])}
    pathNum = 1
    pathsDict = {}
    for startkey, startvalue in startNodeDict.items():
        attempts = 1
        print(attempts)
        finaledgeList = []
        for endkey, endvalue in endNodeDict.items():
            try:
                shortestPath = nx.shortest_path(G, startvalue, endvalue[0], weight='weight')  # (G, source, target)
                if any(point in list(values[0] for values in endNodeDict.values()) for point in shortestPath[1:-1]) or\
                        any(point in list(values for values in barrierNodeDict.values()) for point in shortestPath):
                    print('bad line')
                    pass
                else:
                    pathShape = LineString(shortestPath)  # remove if linestring doesn't work
                    pathsDict[pathNum] = pathShape, pathShape.length, startkey, endkey
                    pathNum += 1
                    # Get network edges that touch the path
                    edgeList = [(u, v) for u, v in G.edges(shortestPath)]
                    print(edgeList)
                    # Remove network edges that are underneath the path
                    edges_to_remove = [edge for edge in edgeList if any(
                        edge[0] in pathShape.coords and edge[1] in pathShape.coords for pathShape, _, _, _ in
                        pathsDict.values())]
                    for edge in edges_to_remove:
                        edgeList.remove(edge)
                    print(edgeList)
                    # Remove network edges that touch an end feature, leaving only open available routes
                    for u, v in edgeList[:]:
                        if endvalue[0] in [u, v]:
                            edgeList.remove((u, v))
                        else:
                            edge_line = LineString([u, v])
                    for each in edgeList:
                        if each not in finaledgeList:
                            finaledgeList.append(each)
                    print(finaledgeList)
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
                print(len(finaledgeList))
                # if finaledge list is at 1
                if len(finaledgeList) >= 1:
                    attempts += 1
                    print(attempts)
                if len(finaledgeList) == 0:
                    print("fully isolated!@!@!@!")
                    arcpy.AddMessage("fully isolated!@!@!@!")
                    break
                if attempts == 1000:
                    print("could not fully isolate. check route layer for open ends.")
                    arcpy.AddWarning("could not fully isolate. check route layer for open ends.")
                    break
            except nx.exception.NetworkXNoPath or nx.exception.NodeNotFound:
                pass

    return pathsDict



def flow_rate(v, diameterFT):

    # Q = v * A
    A = 0.785 * diameterFT**2
    Q = v * A

    # convert to gal/m
    gal = Q * 7.48

    return Q, gal

def fanning_friction_factor(diameter_meters):

    # hf = f(L / Rh) * (v^2/2g)

    Re = 4000  # Reynolds number for turbulent flow
    Re /= 1000  # Convert to thousands

    # Moody chart correlation (for turbulent flow)
    if Re < 2000:
        return 16 / Re
    else:
        f_guess = 0.02  # Initial guess for f
        while True:
            f = 1 / (-2 * math.log10((e / (3.7 * pipeDiameter)) + (2.51 / (Re * math.sqrt(f_guess))))**2)
            if abs(f - f_guess) < 1e-6:
                return f


def darcy_weisbach_headloss(f, L,  D, v, g):
    # final output, h, is in meters

    # h = f (L/D) * (v^2/2g)

    pass
    # return h


def darcy_weisbach_velocity(startingPressure, D, L, f):

    # Calculate velocity using the Darcy-Weisbach equation
    v = math.sqrt((2 * delta_P) / (p * f * length_meters))

    return v


def water_Pressure():
    # P = (p * g * h) + (4 * Q/pi * d^2)
    pass

# Equation Varibles:
pipeMaterial = 'PVC'  # example
startingPressure = 60  # example
pipeDiameter = 4  # example
pipeLength = 1000   # example

# (All values are converted to Metric)
roughnessDict = {'PVC': 0.0015, 'DIP': 0.045, 'HDPE': 0.0025, 'AC': 0.3, 'CAS': 0.26}
psi_to_pa = 6894.76  # Conversion factor from psi to Pascals
p = 1000  # Density of water in kg/m³
g = 9.81  # Acceleration due to gravity in m/s²
e = roughnessDict[pipeMaterial]  # Internal Roughness of Pipe
diameter_meters = pipeDiameter * 0.0254
diameterFT = pipeDiameter / 12
length_meters = pipeLength * 0.3048
delta_P = startingPressure * psi_to_pa

f = fanning_friction_factor(diameter_meters)

print(f"Darcy friction factor: {f}")
v = darcy_weisbach_velocity(startingPressure, diameter_meters, length_meters, f)
print(f"Velocity of water: {v} m/s")
Q, gal = flow_rate(v, diameterFT)
print(f"Flow rate of water: {gal} gpm")

# h = darcy_weisbach_headloss(f, length_meters, diameter_meters, v, g)


def createPathLayer(pathsDict):

    # Create a new polyline feature class
    if pathsDict:
        arcpy.AddMessage("Creating Route layer")
        # arcpy.env.addOutputsToMap = True
        pathFC = arcpy.CreateFeatureclass_management(ws, 'Route', 'POLYLINE', spatial_reference=arcpy.Describe(startFeature).SpatialReference)
        arcpy.management.AddFields(
            in_table=pathFC,
            field_description="Start_Feature TEXT Start_Feature 255 # #;End_Feature TEXT End_Feature 255 # #",
            template=None
        )

        allPathdesc = arcpy.Describe(pathFC)
        dirname = os.path.dirname(arcpy.Describe(pathFC).catalogPath)
        with arcpy.da.Editor(dirname, multiuser_mode=allPathdesc.isVersioned):
            with arcpy.da.InsertCursor(pathFC, ['Start_Feature','End_Feature','SHAPE@']) as cursor:
                for item in pathsDict.values():  # newPathList
                    pathLine = item[0]
                    # turn shapely LineString into an arcgis geometry
                    arcgis_geom = arcpy.FromWKT(pathLine.wkt)
                    # then add it to feature class
                    cursor.insertRow([item[2], item[3], arcgis_geom])

        # if 'Route' not in currentMap.listLayers():
        currentMap.addDataFromPath(pathFC)
    else:
        arcpy.AddMessage("No paths created.")


if __name__ == '__main__':

    import arcpy
    import os
    import sys
    import math
    import networkx as nx
    import pickle
    from shapely.geometry import Point, LineString
    from shapely import STRtree
    from shapely.ops import split
    # from networkx.exception import NetworkXNoPath, NodeNotFound

    arcpy.SetLogMetadata(False)
    arcpy.SetLogHistory(False)
    arcpy.env.overwriteOutput = True
    # sys.tracebacklimit = 0

    RouteType = arcpy.GetParameter(0)

    startFeature = arcpy.GetParameterAsText(1)
    startIDField = arcpy.GetParameterAsText(2)

    endFeature = arcpy.GetParameterAsText(3)
    endIDField = arcpy.GetParameterAsText(4)

    startFeatureFilter = None
    endFeatureFilter = None
    if startFeature == endFeature:
        startFeatureFilter = arcpy.GetParameter(5)
        endFeatureFilter = arcpy.GetParameter(6)

    barrierFeature = arcpy.GetParameterAsText(7)
    barrierFeatureFilter = None
    if barrierFeature:
        barrierFeatureFilter = arcpy.GetParameter(6)

    aprx = arcpy.mp.ArcGISProject('current')
    currentMap = aprx.activeMap
    arcpy.env.workspace = aprx.defaultGeodatabase
    ws = arcpy.env.workspace

    # Load the graph from the file
    with open("graph.pkl", 'rb') as f:
        G = pickle.load(f)

    # try:
    startNodeDict, endNodeDict, graphNodeList, graphNodeTree, barrierNodeDict = getLocations()
    # except IndexError:
    #     # sys.exit()
    #     # sys.tracebacklimit = 1000
    #     raise ValueError("The Start Feature was not found on a vertex. Please snap the feature to a vertex on a related line.")

    if startFeature and endFeature:
        if RouteType == 'Isolation':
            pathsDict = isolationPath(startNodeDict, endNodeDict, graphNodeList, graphNodeTree, barrierNodeDict)
        if RouteType == 'Nearest End Feature':
            pathsDict = createShortestPath(startNodeDict, endNodeDict, graphNodeList, graphNodeTree, barrierNodeDict)

        createPathLayer(pathsDict)

    del aprx, currentMap, G

