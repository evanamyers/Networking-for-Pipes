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

lineNetworkFeatures = arcpy.GetParameterAsText(1).replace("'","")
listLineFC = lineNetworkFeatures.split(";") if ";" in lineNetworkFeatures else [lineNetworkFeatures]

junctions = arcpy.GetParameterAsText(3).replace("'","")  # Should be a Fitting feature
junctionsFC = junctions.split(";") if ";" in junctions else [junctions]

arcpy.AddMessage(str(listLineFC))

lineFilter = arcpy.GetParameterAsText(2)
junctionFilter = arcpy.GetParameterAsText(4)


def createGraph():

    lineNetworkDict = {}
    lineCoordList = []
    lineShapeList = []
    for l, each in enumerate(listLineFC):
        letter = chr(ord('a') + l)
        i = 0
        with arcpy.da.SearchCursor(each, ["GlobalID", "SHAPE@", "SHAPE@LENGTH"], where_clause=lineFilter) as cur:
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
                lineNetworkDict[row[0]] = (nodes, lineXY)

    # Add nodes and edges
    for node_id, (labels, coordinates) in lineNetworkDict.items():
        for label, coord in zip(labels, coordinates):
            G.add_node(coord, id=node_id, label=label)

        for i in range(len(coordinates) - 1):
            x1, y1 = coordinates[i]
            x2, y2 = coordinates[i + 1]
            # add weights to each edge based on length
            G.add_edge(coordinates[i], coordinates[i + 1], weight=math.sqrt((x2 - x1)**2 + (y2 - y1)**2))

    if G:

        edgesCoords = [edge for edge in G.edges()]

        if edgesCoords:
            output_fcEdges = arcpy.CreateFeatureclass_management(arcpy.env.workspace, 'network_edges', 'POLYLINE',
                                                                 spatial_reference=arcpy.Describe(
                                                                     listLineFC[0]).SpatialReference)
            dirname = os.path.dirname(arcpy.Describe(output_fcEdges).catalogPath)
            with arcpy.da.Editor(dirname, multiuser_mode=True):
                with arcpy.da.InsertCursor(output_fcEdges, ['SHAPE@']) as cursor:
                    for coords in edgesCoords:
                        line = LineString(coords)
                        polyline = arcpy.FromWKT(line.wkt)  # Unpack the tuple
                        cursor.insertRow([polyline])

            currentMap.addDataFromPath(edgesCoords)
        arcpy.AddMessage('Created Graph')

    return G


def fixBrokenConnections():
    # try:
    print("Graph is not connected.")

    # Get the connected components
    connected_components = list(nx.connected_components(G))
    largest_component = max(connected_components, key=len)  # this lists just the nodes
    largest_component_edges = G.subgraph(largest_component).edges()  # this lists the edges

    # Find nodes that are not part of the largest connected component
    isolated_nodesCoords = [node for node in G.nodes() if node not in largest_component]
    isolated_edgesCoords = [edge for edge in G.edges() if edge not in largest_component_edges]

    # print("Isolated nodes:", isolated_nodes)

    # print(f"There are: {len(isolated_nodes)} nodes that are isolated")
    # arcpy.AddMessage(f"There are: {len(isolated_nodes)} nodes that are isolated")
    # print(f"There are: {len(isolated_edgesCoords)} edges that are isolated")
    # arcpy.AddMessage(f"There are: {len(isolated_edgesCoords)} edges that are isolated")

    if isolated_nodesCoords:
        output_fcNodes = arcpy.CreateFeatureclass_management(arcpy.env.workspace, 'isolated_nodes', 'POINT',
                                                             spatial_reference=arcpy.Describe(
                                                                 listLineFC[0]).SpatialReference)
        dirname = os.path.dirname(arcpy.Describe(output_fcNodes).catalogPath)
        with arcpy.da.Editor(dirname, multiuser_mode=True):
            with arcpy.da.InsertCursor(output_fcNodes, ['SHAPE@']) as cursor:
                for coords in isolated_nodesCoords:
                    point = arcpy.Point(*coords)  # Unpack the tuple
                    cursor.insertRow([point])

        currentMap.addDataFromPath(output_fcNodes)


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
    if isConnected == False:
        arcpy.AddMessage('Some parts of Network is Disconnected.  Creating an Isolated Lines layer to find areas not '
                         'connected to larger network.')
        fixBrokenConnections()

else:
    arcpy.AddMessage("Creating Single Direction Network")
    print("Creating Single Direction Network")
    G = nx.DiGraph()  # Di meaning 'Directed' (only one direction)  # used to be MultiDiGraph
    createGraph()

    is_strongly_connected = nx.is_strongly_connected(G)
    is_weakly_connected = nx.is_weakly_connected(G)
    arcpy.AddMessage(f"Strongly connected: {is_strongly_connected}")
    arcpy.AddMessage(f"Weakly connected: {is_weakly_connected}")
    print(f"Strongly connected: {is_strongly_connected}")
    print(f"Weakly connected: {is_weakly_connected}")


print(G.edges())
print(G.nodes())

# Save the graph to a file
with open('graph.pkl', 'wb') as f:
    pickle.dump(G, f)


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



    f = fanning_friction_factor(diameter_meters)

    print(f"Darcy friction factor: {f}")
    v = darcy_weisbach_velocity(startingPressure, diameter_meters, length_meters, f)
    print(f"Velocity of water: {v} m/s")
    Q, gal = flow_rate(v, diameterFT)
    print(f"Flow rate of water: {gal} gpm")
    # h = darcy_weisbach_headloss(f, length_meters, diameter_meters, v, g)
