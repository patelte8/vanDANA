fsi_interpolation_code = """

/*
<%
from dolfin.jit.jit import dolfin_pc
setup_pybind11(cfg)
cfg['include_dirs'] = dolfin_pc['include_dirs']
cfg['library_dirs'] = dolfin_pc['library_dirs']
cfg['compiler_args']  = ['-std=c++11', '-DHAS_MPI']
%>
*/

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <time.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/common/RangedIndexSet.h>
#include <dolfin/geometry/BoundingBoxTree.h>
#include <dolfin/mesh/Edge.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/CellType.h>
#include <dolfin/mesh/MeshEntityIterator.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/geometry/Point.h>
#include <dolfin/mesh/MeshEntity.h>

using namespace dolfin;
namespace py = pybind11;

std::unordered_map<std::size_t, std::size_t> dof_component_map_f1, dof_component_map_f2, dof_component_map_s1, dof_component_map_s2;

/// Comparison operator for hashing coordinates. Note that two
// coordinates are considered equal if equal to within specified
// tolerance.
struct lt_coordinate
{
  lt_coordinate(double tolerance) : TOL(tolerance) {}

  bool operator() (const std::vector<double>& x,
                   const std::vector<double>& y) const
  {
    std::size_t n = std::max(x.size(), y.size());
    for (std::size_t i = 0; i < n; ++i)
    {
      double xx = 0.0;
      double yy = 0.0;
      if (i < x.size())
        xx = x[i];
      if (i < y.size())
        yy = y[i];

      if (xx < (yy - TOL))
        return true;
      else if (xx > (yy + TOL))
        return false;
    }
    return false;
  }

  // Tolerance
  const double TOL;
};

void extract_dof_component_map(std::unordered_map<std::size_t,
                               std::size_t>& dof_component_map,
                               const FunctionSpace& V,
                               int* component)
{
  // Extract sub dofmaps recursively and store dof to component map
  if (V.element()->num_sub_elements() == 0)
  {
    std::unordered_map<std::size_t, std::size_t> collapsed_map;
    std::shared_ptr<GenericDofMap> dummy
      = V.dofmap()->collapse(collapsed_map, *V.mesh());
    (*component)++;
    for (const auto &map_it : collapsed_map)
      dof_component_map[map_it.second] = (*component);
  }
  else
  {
    for (std::size_t i = 0; i < V.element()->num_sub_elements(); ++i)
    {
      const std::vector<std::size_t> comp = {i};
      std::shared_ptr<FunctionSpace> Vs = V.extract_sub_space(comp);
      extract_dof_component_map(dof_component_map, *Vs, component);
    }
  }}

void extract_dof_component_map_user(const FunctionSpace& Z, std::string flag)
{
  int component = -1;
  const FiniteElement& element = *Z.element();

  if (flag == "F") {
    if (element.value_rank() == 1)
      extract_dof_component_map(dof_component_map_f1, Z, &component);
    if (element.value_rank() == 0)
      extract_dof_component_map(dof_component_map_f2, Z, &component);
  }

  if (flag == "S") {
    if (element.value_rank() == 1)
      extract_dof_component_map(dof_component_map_s1, Z, &component);
    if (element.value_rank() == 0)
      extract_dof_component_map(dof_component_map_s2, Z, &component);
  }
}

bool in_bounding_box(const std::vector<double>& point,
                     const std::vector<double>& bounding_box,
                     const double tol)
{
  // Return false if bounding box is empty
  if (bounding_box.empty())
    return false;

  const std::size_t gdim = point.size();
  dolfin_assert(bounding_box.size() == 2*gdim);
  for (std::size_t i = 0; i < gdim; ++i)
  {
    if (!(point[i] >= (bounding_box[i] - tol)
          && point[i] <= (bounding_box[gdim + i] + tol)))
    {
      return false;
    }
  }
  return true;}


std::vector<double> create_mini_bounding_box(const std::vector<double>& point, double h, const double tol)
{
  const std::size_t gdim = point.size();
  std::vector<double> m_m(2*gdim);  

  for (std::size_t i = 0; i < gdim; ++i)
  {
        m_m[i]        = point[i] - (3*h + tol);
        m_m[gdim + i] = point[i] + (3*h + tol);
  }

  return m_m;
}


std::vector<double> create_bounding_box(const MPI_Comm mpi_comm, const Mesh& mesh0, const Mesh& mesh1, double& h, const std::size_t gdim, std::string meshid)
{
  std::vector<double> coordinates;

  if (meshid == "F")
    coordinates = mesh0.coordinates();
  if (meshid == "S")
    coordinates = mesh1.coordinates();
  
  std::size_t num_processes = dolfin::MPI::size(mpi_comm);
  std::vector<double> x_min_max(2*gdim);
  
  for (std::size_t i = 0; i < gdim; ++i)
  {
      x_min_max[i]        = *(coordinates.begin() + i);
      x_min_max[gdim + i] = *(coordinates.begin() + i);
  }

  for (std::size_t i = 0; i < gdim; ++i)
  {
      for (auto it = coordinates.begin() + i; it < coordinates.end(); it += gdim)
      {
      x_min_max[i]        = std::min(x_min_max[i], *it);
      x_min_max[gdim + i] = std::max(x_min_max[gdim + i], *it);
      }
  }

  // Communicate bounding boxes
  std::vector<std::vector<double>> bounding_boxes;
  dolfin::MPI::all_gather(mpi_comm, x_min_max, bounding_boxes);
  
  // Create overall bounding box for solid mesh across processes
  std::vector<double> x_1(2*gdim);
  for (std::size_t i = 0; i < gdim; ++i)
  {
    for (auto it = bounding_boxes[0].begin() + i; it < bounding_boxes[0].end(); it += gdim)
    {
        x_1[i]        = *it;
        x_1[gdim + i] = *it;
    }
  }

  for(std::size_t p = 0; p < num_processes; p++)
  {
    for (std::size_t i = 0; i < gdim; ++i)
    {
        for (auto it = bounding_boxes[p].begin() + i; it < bounding_boxes[p].end(); it += gdim)
        {
            x_1[i]        = std::min(x_1[i], *it);
            x_1[gdim + i] = std::max(x_1[gdim + i], *it);
        }
    }
  }

  // Calculate eulerian meshsize
  std::vector<double> x(gdim);
  double h_local;
  double minh = 10000;
  double maxh = -10000;
  
  if (meshid == "F")
  { 
      for (EdgeIterator edge(mesh1); !edge.end(); ++edge) 
      {
          for (VertexIterator vert(*edge); !vert.end(); ++vert)
          {   
              for(std::size_t i = 0; i < gdim; i++)
                  x[i]=vert->x(i); 
  
              if(in_bounding_box(x, bounding_boxes[dolfin::MPI::rank(mpi_comm)], 1e-12))
              { 
                  if(edge->length() < minh) {minh = edge->length();};
                  if(edge->length() > maxh) {maxh = edge->length();};
              }
          }
      }
  }
  
  else if (meshid == "S")
  {
      for (EdgeIterator edge(mesh0); !edge.end(); ++edge)
      {
          for (VertexIterator vert(*edge); !vert.end(); ++vert)
          {    
              for(std::size_t i = 0; i < gdim; i++)
                  x[i]=vert->x(i);

              if(in_bounding_box(x, bounding_boxes[dolfin::MPI::rank(mpi_comm)], 1e-12))
              { 
                  if(edge->length() < minh) {minh = edge->length();};
                  if(edge->length() > maxh) {maxh = edge->length();};
              }
          }   
      }
  }
  
  h_local = minh;
  h = dolfin::MPI::min(mpi_comm, h_local);

  // Extending bounding box of mesh by 3*h to determine points for MPI
  for (std::size_t i = 0; i < gdim; ++i)
  {
      x_1[i]        -= 3.0*h;
      x_1[gdim + i] += 3.0*h;
  }

  return x_1;
}


double phi1(double frac)
{
  double value = 0.0;

  if(fabs(frac) <= 1.0)
      value = 1 - fabs(frac);
  
  return value;
}

double phi1_smoothed(double frac)
{
  double value = 0.0;

  if(fabs(frac) <= 0.5)
      value = (3.0/4) - (frac*frac);
  else if(fabs(frac) >= 0.5 && fabs(frac) <= 1.5)
      value = (9.0/8) - ((3.0/2)*fabs(frac)) + ((frac*frac)/2.0);

  return value;
}

double phi2(double frac)
{
  double value = 0.0;

  if(fabs(frac) <= 2.0)
      value = 0.25*(1+cos(0.5*3.14159265359*frac));

  return value;
}

double phi2_smoothed(double frac)
{
  double value = 0.0;

  if(fabs(frac) <= 1.5)
      value = (1/(4*3.14159265359))*(3.14159265359 + (2*sin((3.14159265359/4)*((2*frac)+1))) - (2*sin((3.14159265359/4)*((2*frac)-1))));
  else if(fabs(frac) >= 1.5 && fabs(frac) <= 2.5)
      value = (-1/(8*3.14159265359))*(((-1*5)*3.14159265359) + (2*3.14159265359*fabs(frac)) + (4*sin((3.14159265359/4)*((2*fabs(frac)-1)))));
  
  return value;
}

double phi3(double frac)
{
  double value = 0.0;

  if(fabs(frac) <= 0.5)
      value = (1.0/3)*(1.0+std::sqrt((-3*frac*frac)+1.0));
  else if(fabs(frac) >= 0.5 && fabs(frac) <= 1.5)
      value = (1.0/6)*(5.0-(3*fabs(frac))-std::sqrt(((-1*3.0)*(1-fabs(frac))*(1-fabs(frac)))+1.0));

  return value;
}

double phi3_smoothed(double frac)
{
  double value = 0.0;

  if(fabs(frac) <= 1.0)
      value = (17.0/48) + ((1.732050808*3.14159265359)/108) + (fabs(frac)/4) - ((frac*frac)/4) + 
  			(((1.0-(2*(fabs(frac))))/16)*(std::sqrt((-1*12*frac*frac) + (12*fabs(frac)) + 1.0))) 
  			- ((1.732050808/12)*(asin((1.732050808/2)*((2*fabs(frac))-1))));
  else if(fabs(frac) >= 1.0 && fabs(frac) <= 2.0)
      value = (55.0/48) - ((1.732050808*3.14159265359)/108) - ((13*fabs(frac))/12) + ((frac*frac)/4) + 
  			((((2*(fabs(frac)))-3.0)/48)*(std::sqrt((-1*12*frac*frac) + (36*fabs(frac)) - 23.0)))
  			+ ((1.732050808/36)*(asin((1.732050808/2)*((2*fabs(frac))-3))));

  return value;
}

double phi4(double frac)
{
  double value = 0.0;

  if(fabs(frac) <= 1.0)
      value = (0.125)*(3.0-(2*fabs(frac))+std::sqrt(1.0+(4*fabs(frac))-(4*frac*frac)));
  else if(fabs(frac) >= 1.0 && fabs(frac) <= 2.0)
      value = (0.125)*(5.0-(2*fabs(frac))-std::sqrt((-1*7.0)+(12*fabs(frac))-(4*frac*frac)));
  
  return value;
}

double phi4_smoothed(double frac)
{
  double value = 0.0;;

  if(fabs(frac) <= 0.5)
      value = (3.0/8) + (3.14159265359/32) - ((fabs(frac)*fabs(frac))/4);
  else if(fabs(frac) >= 0.5 && fabs(frac) <= 1.5)
      value = (1.0/4) + (((1.0-fabs(frac))/8)*(std::sqrt((-1*4.0*frac*frac) + (8.0*fabs(frac)) - 2.0)))
  			- ((1.0/8)*(asin(1.414213562*(fabs(frac)-1))));
  else if(fabs(frac) >= 1.5 && fabs(frac) <= 2.5)
      value = (17.0/16) - (3.14159265359/64) - ((3*fabs(frac))/4) + ((frac*frac)/8) + 
  			(((fabs(frac)-2.0)/16)*(std::sqrt((-1*4.0*frac*frac) + (16*fabs(frac)) - 14.0)))
  			+ ((1.0/16)*(asin(1.414213562*(fabs(frac)-2))));
  
  return value;
}

double phi5(double frac)
{
  double value = 0.0;

  if(fabs(frac) <= 2.0)
      value = 0.5 - (0.25*fabs(frac));
  
  return value;
}

const static std::unordered_map<std::string,int> string_to_case{
   {"phi1",1},
   {"phi2",2},
   {"phi3",3},
   {"phi4",4},
   {"phi5",5},
   {"phi1_smoothed",6},
   {"phi2_smoothed",7},
   {"phi3_smoothed",8},
   {"phi4_smoothed",9}
};

double delta(std::size_t gdim, double& h, Array<double>& x, std::vector<double>& node_coord, std::string& abc)
{
  double wax = 1;  

  switch (string_to_case.at(abc)) {
    
    case 1:

        for(std::size_t i = 0; i < gdim; i++)      
          {wax *= phi1((x[i] - node_coord[i])/h)/h;}
        break;
 
    case 2:

        for(std::size_t i = 0; i < gdim; i++)      
          {wax *= phi2((x[i] - node_coord[i])/h)/h;}
        break;

    case 3:

        for(std::size_t i = 0; i < gdim; i++)      
          {wax *= phi3((x[i] - node_coord[i])/h)/h;}
        break;
          
    case 4:

        for(std::size_t i = 0; i < gdim; i++)      
          {wax *= phi4((x[i] - node_coord[i])/h)/h;}
        break;

    case 5:

        for(std::size_t i = 0; i < gdim; i++)      
          {wax *= phi5((x[i] - node_coord[i])/h)/h;}
        break;    
          
    case 6:

        for(std::size_t i = 0; i < gdim; i++)      
          {wax *= phi1_smoothed((x[i] - node_coord[i])/h)/h;}
        break;
          
    case 7:

        for(std::size_t i = 0; i < gdim; i++)      
          {wax *= phi2_smoothed((x[i] - node_coord[i])/h)/h;}
        break;

    case 8:

        for(std::size_t i = 0; i < gdim; i++)      
          {wax *= phi3_smoothed((x[i] - node_coord[i])/h)/h;}
        break;
          
    case 9:

        for(std::size_t i = 0; i < gdim; i++)      
          {wax *= phi4_smoothed((x[i] - node_coord[i])/h)/h;}
        break;                        
  }

  return wax;
}

double fraction(const CellType::Type celltype)
{
  double frac = 0.0;
  // double std_vol;   

  if(celltype == CellType::Type::triangle)
      frac = 1.0/3;
      //std_vol = 0.5;    

  else if(celltype == CellType::Type::tetrahedron)
      frac = 1.0/4;
      //std_vol = 1.0/6; 

  else if(celltype == CellType::Type::hexahedron)
      frac = 1.0/6;
      //std_vol = 1.0;

  else if(celltype == CellType::Type::quadrilateral)
      frac = 1.0/4;
      //std_vol = 1.0; 

  else
      {std::cout << "Mesh celltype not recognised";}

  return frac;
}

std::map<std::vector<double>, std::vector<std::size_t>, lt_coordinate>
tabulate_coordinates_to_dofs(const FunctionSpace& V, std::vector<double>& x_y)
{
    std::map<std::vector<double>, std::vector<std::size_t>, lt_coordinate>
    coords_to_dofs(lt_coordinate(1.0e-12));

  // Extract mesh, dofmap and element
  const GenericDofMap& dofmap = *V.dofmap();
  const FiniteElement& element = *V.element();
  const Mesh& mesh = *V.mesh();

  // Geometric dimension
  const std::size_t gdim = mesh.geometry().dim();

  // Loop over cells and tabulate dofs
  boost::multi_array<double, 2> coordinates;
  std::vector<double> coordinate_dofs;
  std::vector<double> coors(gdim);

  // Speed up the computations by only visiting (most) dofs once
  const std::size_t local_size = dofmap.ownership_range().second
    - dofmap.ownership_range().first;
  RangedIndexSet already_visited(std::make_pair(0, local_size));

  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    // Update UFC cell
    cell->get_coordinate_dofs(coordinate_dofs);

    // Get local-to-global map
    auto dofs = dofmap.cell_dofs(cell->index());

    // Tabulate dof coordinates on cell
    element.tabulate_dof_coordinates(coordinates, coordinate_dofs,
                                     *cell);

    // Map dofs into coords_to_dofs
    for (Eigen::Index i = 0; i < dofs.size(); ++i)
    {
      const std::size_t dof = dofs[i];
      if (dof < local_size)
      {
        // Skip already checked dofs
        if (!already_visited.insert(dof))
          continue;

        // Put coordinates in coors
        std::copy(coordinates[i].begin(), coordinates[i].end(),
                  coors.begin());

        // Add dof to list at this coord
        if(in_bounding_box(coors, x_y, 1e-12))
        {
          const auto ins = coords_to_dofs.insert
            (std::make_pair(coors, std::vector<std::size_t>{dof}));
          if (!ins.second)
            ins.first->second.push_back(dof);
        }
      }
    }
  }
  return coords_to_dofs;
}

std::map<std::vector<double>, std::pair<std::vector<double>, double>>
tabulate_coordinates_to_values_and_nodal_volumes(const Function& u0, std::vector<double>& x_x)
{
	std::map<std::vector<double>, std::pair<std::vector<double>, double>>
  	coords_to_values_and_nodal_volumes;        

  // Extract mesh
  const Mesh& mesh0 = *u0.function_space()->mesh();
  const CellType::Type celltype = mesh0.type().cell_type();
  const std::size_t gdim0 = mesh0.geometry().dim();

  // Create arrays used to evaluate points on mesh0
  double nodal_volume;
  std::vector<double> node_coord(gdim0);
	std::vector<double> node_val(u0.value_size());
	Array<double> _node_coord(gdim0, node_coord.data());
	Array<double> _node_val(u0.value_size(), node_val.data());

  for (VertexIterator vert(mesh0); !vert.end(); ++vert)    
	{
 	  for(std::size_t i = 0; i < gdim0; i++)
      node_coord[i] = vert->x(i);

    u0.eval(_node_val, _node_coord);      // Remember this is a known bug here.

    if(in_bounding_box(node_coord, x_x, 1e-12))
    {
  	nodal_volume = 0;
  	for (CellIterator cell(*vert); !cell.end(); ++cell)
  	    nodal_volume += (*cell).volume();

  	nodal_volume *= fraction(celltype);

  	coords_to_values_and_nodal_volumes.insert
  	  (std::make_pair(node_coord, std::make_pair(node_val, nodal_volume)));	
    }
  }

return coords_to_values_and_nodal_volumes;
}


void interpolate_delta(const Function& u0, Function& u, std::string abc, std::string meshid) 
{
  // Interpolate from Function u0 to Function u.
  // This mesh of u0 may be different from that of u
  //
  // The algorithm is briefly
  //
  //   1) Create a map from all different coordinates of u dofs to
  //      the dofs living on that coordinate. This is done such that
  //      one only need to visit (and distribute) each interpolation
  //      point once.
  //   2) Create a map from dof to component index in Mixed Space.
  //   3) Create bounding boxes for the partitioned mesh of u0 and
  //      distribute to all processors.
  //   4) Using bounding boxes, compute the processes that *may* own
  //      the dofs of u.
  //   5) Distribute interpolation points to potential owners who
  //      subsequently tries to evaluate u0. If successful, return
  //      values of u0 to owner.

  // Get function spaces of Functions interpolating to/from
  dolfin_assert(u0.function_space());
  dolfin_assert( u.function_space());
  const FunctionSpace& V0 = *u0.function_space();
  const FunctionSpace& V1 =  *u.function_space();

  // Get element interpolating to
  dolfin_assert(V1.element());
  const FiniteElement& element = *V1.element();

  // Check that function ranks match
  if (element.value_rank() != u0.value_rank())
  {
      dolfin_error("LagrangeInterpolator.cpp",
                  "interpolate Function into function space",
                  "Rank of Function (%d) does not match rank of function space (%d)",
                  u0.value_rank(), element.value_rank());
  }

  // Check that function dims match
  for (std::size_t i = 0; i < element.value_rank(); ++i)
  {
      if (element.value_dimension(i) != u0.value_dimension(i))
      {
      dolfin_error("LagrangeInterpolator.cpp",
                  "interpolate Function into function space",
                  "Dimension %d of Function (%d) does not match dimension %d of function space (%d)",
                  i, u0.value_dimension(i), i, element.value_dimension(i));
      }
  }

  // Get mesh and dimension of FunctionSpace interpolating to/from
  dolfin_assert(V0.mesh());
  dolfin_assert(V1.mesh());
  const Mesh& mesh0 = *V0.mesh();
  const Mesh& mesh1 = *V1.mesh();
  const std::size_t gdim0 = mesh0.geometry().dim();

  // Get communicator
  const MPI_Comm mpi_comm = V1.mesh()->mpi_comm();

  // Get number of MPI processes
  std::size_t num_processes = dolfin::MPI::size(mpi_comm);

  // Create arrays used to evaluate points on mesh1
  std::vector<double> x(gdim0);
  std::vector<double> values(u0.value_size());
  Array<double> _x(gdim0, x.data());
  Array<double> _values(u0.value_size(), values.data());

  // Create arrays used to evaluate points on mesh0
  std::vector<double> node_coord(gdim0);
  std::vector<double> node_val(u0.value_size());
  Array<double> _node_coord(gdim0, node_coord.data());
  Array<double> _node_val(u0.value_size(), node_val.data());

  // Create bounding box for solid mesh 
  std::vector<double> x_1(2*gdim0); double h = 0.0;
  x_1 = create_bounding_box(mpi_comm, mesh0, mesh1, h, gdim0, meshid);

  // Create vector to hold all local values of u
  std::vector<double> local_u_vector(u.vector()->local_size());

  // Create map from coordinates to dofs sharing that coordinate
  std::map<std::vector<double>, std::vector<std::size_t>, lt_coordinate>
      coords_to_dofs = tabulate_coordinates_to_dofs(V1, x_1);
  
  // Map from coordinates to respective values and nodal volumes
  double nodal_volume;
  std::vector<double> m_m(2*gdim0); 
  std::map<std::vector<double>, std::pair<std::vector<double>, double>> 
      coords_to_values_and_nodal_volumes = tabulate_coordinates_to_values_and_nodal_volumes(u0, x_1);    

  // Search this process first for all coordinates in u local mesh
  std::vector<std::vector<double>> all_points(num_processes);
	for (std::size_t p = 0; p < num_processes; p++)
  { 
		for (const auto &map_it : coords_to_dofs)
		{
		    // Place interpolation point in x
		    std::copy(map_it.first.begin(), map_it.first.end(), x.begin());
		  
		   	if (in_bounding_box(x, x_1, 1e-12))
				all_points[p].insert(all_points[p].end(), x.begin(), x.end());		    
		}
	}

  // Communicate all potential points
  std::vector<std::vector<double>> potential_points_recv;
  dolfin::MPI::all_to_all(mpi_comm, all_points, potential_points_recv);

  // Now try to eval u0 for the received points
  std::vector<std::vector<double>> coefficients_found(num_processes);
  std::vector<std::vector<double>> points_found(num_processes);
  for (std::size_t p = 0; p < num_processes; ++p)
  {
      std::vector<double>& points = potential_points_recv[p];
      for (std::size_t j = 0; j < points.size()/gdim0; ++j)
      {
        std::copy(points.begin() + j*gdim0, points.begin() + (j + 1)*gdim0,
                x.begin());
	    
		    // ................ Delta function interpolation ..................
		    m_m = create_mini_bounding_box(x, h, 1e-12);

		    // Initializing values as zero
		    for(std::size_t i = 0; i < gdim0; i++)
		        _values[i] = 0;
		    
		    for (const auto &map : coords_to_values_and_nodal_volumes)
			  {
  		    	std::copy(map.first.begin(), map.first.end(), node_coord.begin());

  		      if (in_bounding_box(node_coord, m_m, 1e-12))
  		      {	
              auto pr = map.second;
  		      	std::copy(pr.first.begin(), pr.first.end(), node_val.begin());
              nodal_volume = pr.second;

  		       	for(std::size_t i = 0; i < gdim0; i++)
  		         	_values[i] += _node_val[i]*delta(gdim0, h, _x, node_coord, abc)*nodal_volume;
  		      }
		    }
		    // ................. Delta function interpolation ..................

        coefficients_found[p].insert(coefficients_found[p].end(),
	                                        values.begin(), values.end());
        points_found[p].insert(points_found[p].end(), x.begin(), x.end());
      }
  }

  // Send back the found coefficients and points
  std::vector<std::vector<double>> coefficients_recv;
  std::vector<std::vector<double>> points_recv;
  dolfin::MPI::all_to_all(mpi_comm, coefficients_found, coefficients_recv);
  dolfin::MPI::all_to_all(mpi_comm, points_found, points_recv);
  for (std::size_t p = 0; p < num_processes; ++p)
  {
      // Get the new values and points
      const std::vector<double>& vals = coefficients_recv[p];
      const std::vector<double>& pts = points_recv[p];

      // Move all found coefficients into the local_u_vector
      for (std::size_t j = 0; j < pts.size()/gdim0; ++j)
      {
      std::copy(pts.begin() + j*gdim0, pts.begin() + (j + 1)*gdim0, x.begin());

      // Get the owned dofs sharing x
      const std::vector<std::size_t>& dofs = coords_to_dofs[x];

      // Place result in local_u_vector
      for (const auto &d : dofs)
      {
          dolfin_assert(d <  local_u_vector.size());
          if (meshid == "S")
          {   
            if (element.value_rank() == 0)
              {local_u_vector[d] += vals[j*u0.value_size() + dof_component_map_s2[d]];}
            if (element.value_rank() == 1)  
              {local_u_vector[d] += vals[j*u0.value_size() + dof_component_map_s1[d]];}
          }
          else if (meshid == "F")
          {
            if (element.value_rank() == 0)
              {local_u_vector[d] += vals[j*u0.value_size() + dof_component_map_f2[d]];}
            if (element.value_rank() == 1)
              {local_u_vector[d] += vals[j*u0.value_size() + dof_component_map_f1[d]];} 
          }
      }            
      }
  }
  
  // Set and finalize
  u.vector()->set_local(local_u_vector);
  u.vector()->apply("insert");
}

PYBIND11_MODULE(SIGNATURE, m)
{
  m.def("interpolate_delta", (void (*)(const Function&, Function&, std::string, std::string))
      &interpolate_delta);
  m.def("interpolate_delta", [](py::object u0, py::object v, py::str abc, py::str meshid){
      auto _u = u0.attr("_cpp_object");
      auto _v = v.attr("_cpp_object").cast<Function*>();
      auto _u0 = _u.cast<const Function*>();
      interpolate_delta(*_u0, *_v, abc.cast<std::string>(), meshid.cast<std::string>());
    });

  m.def("extract_dof_component_map_user", (void (*)(const FunctionSpace &, std::string))
          &extract_dof_component_map_user);
  m.def("extract_dof_component_map_user", [](py::object X, py::str flag){
      auto _X = X.attr("_cpp_object").cast<const FunctionSpace &>();
      extract_dof_component_map_user(_X, flag.cast<std::string>());
  });
}    

""" 
