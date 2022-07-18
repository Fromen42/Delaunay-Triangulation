import numpy as np
from typing import Any
import matplotlib.cm as cm
import matplotlib.pyplot as plt 
from dataclasses import dataclass
from matplotlib import collections  as mc

@dataclass
class Vertex:
    x: float
    y: float
    index: int
    inc_edge:Any

@dataclass
class Half_Edge:
    origin: Vertex
    target: Vertex
    prev_edge: Any
    next_edge: Any
    twin:Any

def plot_dcel_edge(dcel_edge_list,voronoi,vd_edge):
  _, ax = plt.subplots()
  colors = iter(cm.rainbow(np.linspace(0, 1, len(vertexes))))
  # plot delauney
  edge_list_dt = [[np.array([halfedge.origin.x,halfedge.origin.y]),np.array([halfedge.target.x,halfedge.target.y])] for halfedge in dcel_edge_list]
  dt = mc.LineCollection(edge_list_dt, colors=(0.0, 0.75, 0.0, 1),linewidths=2, color=next(colors))
  ax.add_collection(dt)
  ax.autoscale()
  # plot voronoi
  #for vertex in voronoi:
    #plt.scatter(vertex[0],vertex[1])
  #vd = mc.LineCollection(vd_edge, colors=(0.0, 0.75, 0.0, 1),linewidths=2, color=next(colors))
  #ax.add_collection(vd)
  plt.show()


def inCircleTest(a_x: float, a_y: float,
                 b_x: float, b_y: float,
                 c_x: float, c_y: float,
                 d_x: float, d_y: float,):
    matrix = np.array([[a_x, a_y, a_x**2 + a_y**2, 1],
                       [b_x, b_y, b_x**2 + b_y**2, 1],
                       [c_x, c_y, c_x**2 + c_y**2, 1],
                       [d_x, d_y, d_x**2 + d_y**2, 1]]) 
    det = np.linalg.det(matrix)
    if(det > 0):
        return True
    else:
        return False


def define_circle(p1, p2, p3):
    temp = p2[0] * p2[0] + p2[1] * p2[1]
    bc = (p1[0] * p1[0] + p1[1] * p1[1] - temp) / 2
    cd = (temp - p3[0] * p3[0] - p3[1] * p3[1]) / 2
    det = (p1[0] - p2[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p2[1])
    if abs(det) < 1.0e-6:
        return (None, np.inf)
    # Center of circle
    cx = (bc*(p2[1] - p3[1]) - cd*(p1[1] - p2[1])) / det
    cy = ((p1[0] - p2[0]) * cd - (p2[0] - p3[0]) * bc) / det
    radius = np.sqrt((cx - p1[0])**2 + (cy - p1[1])**2)
    return [cx, cy]


def swap_test_dcel(a_index: int,b_index: int,c_index: int):
  global edge_list_dcel
  global vertexes_dcel
  #if np.abs(b_index - c_index) == 1 or np.abs(b_index - c_index) == 6:
  ## The above condition was the issue ##
  if False:
    print("CH boundary, return")
    pass
  else:
    a = vertexes_dcel[a_index]
    b = vertexes_dcel[b_index]
    c = vertexes_dcel[c_index]
    for y in c.inc_edge:
      if y.target == b and y.twin != Any and y.twin.next_edge.next_edge.target == b:
        x = y.twin.next_edge
        d = x.target
        print("a:", a.index,"b:", b.index,"c:", c.index,"d:", d.index)
        #print("edge x origin:",x.origin,"edge y origin:",y.origin)
        if inCircleTest(b.x,b.y,c.x,c.y,d.x,d.y,a.x,a.y):
          print("flip bc:", b.index, c.index, "to add:", a.index, d.index)
          add_edge_dcel_1 = Half_Edge(origin=a,target=d,prev_edge=y.next_edge,next_edge=x.next_edge,twin= Any)
          add_edge_dcel_2 = Half_Edge(origin=d,target=a,prev_edge=x,next_edge=y.prev_edge,twin= add_edge_dcel_1)
          add_edge_dcel_1.twin = add_edge_dcel_2
          # recount incident edge
          c.inc_edge.remove(y)
          b.inc_edge.remove(y.twin)
          a.inc_edge.append(add_edge_dcel_1)
          d.inc_edge.append(add_edge_dcel_2)
          # set dcel for new triangel 2 adc
          x.next_edge = add_edge_dcel_2
          x.prev_edge = y.prev_edge
          y.prev_edge.next_edge = x
          y.prev_edge.prev_edge = add_edge_dcel_2
          # set dcel for new triangel 1 abd
          y.next_edge.next_edge = add_edge_dcel_1
          y.next_edge.prev_edge = y.twin.prev_edge
          y.twin.prev_edge.next_edge = y.next_edge
          y.twin.prev_edge.prev_edge = add_edge_dcel_1
          # remove add edges
          edge_list_dcel.remove(y)
          edge_list_dcel.remove(y.twin)
          edge_list_dcel.append(add_edge_dcel_1)
          edge_list_dcel.append(add_edge_dcel_2)
          print("propagate swap test to adc", a.index, d.index, c.index)
          swap_test_dcel(a.index, d.index, c.index)
          print("propagate swap test to abd", a.index, b.index, d.index)
          swap_test_dcel(a.index, b.index, d.index)
        else:
          print("No need of flip ")


# read example input
# randomly delete vertex and find neighours
vertexes = np.loadtxt("ExampleInput.txt")
#vertexes = np.loadtxt("12.txt")
#vertexes = np.loadtxt("23.txt")
convex_hull = [[vertexes[i],vertexes[(i+1)%len(vertexes)]] for i in range(len(vertexes))]
vertexes_permutated = np.random.permutation(vertexes)
vertexes_permutated_index = [np.where(vertexes==i)[0][0] for i in vertexes_permutated]
neighbours_stack = []
neighbours_stack_id = []
vertexes_index = np.arange(len(vertexes))
for i in vertexes_permutated_index:
  i_index = np.where(vertexes_index==i)[0]
  neighbours_next = vertexes_index[np.asscalar((i_index+1)%len(vertexes_index))]
  neighbours_prev = vertexes_index[np.asscalar(i_index-1)]
  neighbours_stack.append([vertexes[neighbours_prev],vertexes[neighbours_next]])
  neighbours_stack_id.append([neighbours_prev,neighbours_next])
  vertexes_index = np.delete(vertexes_index,i_index)
  if len(vertexes_index) <=3:
    break

# Reverse permuted neighbours
vertexes_permutated_reverse = vertexes_permutated[::-1]
neighbours_stack_reverse = neighbours_stack[::-1]
neighbours_stack_id_reverse = neighbours_stack_id[::-1]
vertexes_permutated_index_reverse = vertexes_permutated_index[::-1]
base_vertex_id = np.array(sorted(vertexes_permutated_index_reverse[:3]))

vertexes_dcel = []
edge_list_dcel = []
# build initial dcel data structure
for vertex_id in range(len(vertexes)):
  vertexes_dcel.append(Vertex(x=vertexes[vertex_id][0],y=vertexes[vertex_id][1],index=vertex_id, inc_edge=[]))
edge_list_dcel = [Half_Edge(origin=vertexes_dcel[vertex1],target = vertexes_dcel[vertrx2],prev_edge=Any, next_edge=Any, twin=Any) for vertex1, vertrx2 in zip(base_vertex_id,np.roll(base_vertex_id,-1))]
for vertex_id in base_vertex_id:
  vertexes_dcel[vertex_id].inc_edge = [edge for edge in edge_list_dcel if edge.origin == vertexes_dcel[vertex_id]]
for edge_id in range(len(edge_list_dcel)):
  edge_list_dcel[edge_id].next_edge = edge_list_dcel[(edge_id+1)%len(edge_list_dcel)]
  edge_list_dcel[edge_id].prev_edge = edge_list_dcel[edge_id-1]

# build delauney triangulation
for id in range(len(vertexes_permutated_index_reverse[3:])):
  #create new dcel items
  add_vertex_dcel = vertexes_dcel[vertexes_permutated_index_reverse[id+3]]
  add_vertex_prev_dcel = vertexes_dcel[neighbours_stack_id_reverse[id][0]]
  add_vertex_next_dcel = vertexes_dcel[neighbours_stack_id_reverse[id][1]]
  add_half_edge_prev_dcel = Half_Edge(origin=add_vertex_prev_dcel,target = add_vertex_dcel, prev_edge=Any, next_edge=Any, twin=Any)
  add_half_edge_next_dcel = Half_Edge(origin=add_vertex_dcel,target = add_vertex_next_dcel, prev_edge=Any, next_edge=Any, twin=Any)
  add_half_edge_oppo_dcel = Half_Edge(origin=add_vertex_next_dcel,target = add_vertex_prev_dcel, prev_edge=Any, next_edge=Any, twin=Any)
  add_edge_list = [add_half_edge_prev_dcel,add_half_edge_next_dcel,add_half_edge_oppo_dcel]
  # build dcel connections
  for i in range(3):
    add_edge_list[i].origin.inc_edge.append(add_edge_list[i])
    add_edge_list[i].prev_edge = add_edge_list[i-1]
    add_edge_list[i].next_edge = add_edge_list[(i+1)%3]
  add_half_edge_oppo_dcel.twin = [x for x in add_half_edge_oppo_dcel.target.inc_edge if x.target == add_half_edge_oppo_dcel.origin][0]
  add_half_edge_oppo_dcel.twin.twin = add_half_edge_oppo_dcel
  # check and flip edges
  edge_list_dcel+= add_edge_list
  swap_test_dcel(add_vertex_dcel.index,add_vertex_prev_dcel.index,add_vertex_next_dcel.index)


# build voronoi diagram
voronoi = []
vd_edge = []

def voronoi_site(dt_edge):
  a = dt_edge.origin
  b = dt_edge.target
  c = dt_edge.next_edge.target
  site = define_circle([a.x,a.y],[b.x,b.y],[c.x,c.y])
  return site

for dt_edge in edge_list_dcel:
  if dt_edge.twin == Any:
    site = voronoi_site(dt_edge)    
    voronoi.append(site)
    vd_edge.append([site, [(dt_edge.origin.x+dt_edge.target.x)/2,(dt_edge.origin.y+dt_edge.target.y)/2]])
  else:
    site1 = voronoi_site(dt_edge)
    site2 = voronoi_site(dt_edge.twin)
    edge = [site1,site2]
    if edge not in vd_edge:
      vd_edge.append(edge)


plot_dcel_edge(edge_list_dcel,voronoi,vd_edge)  