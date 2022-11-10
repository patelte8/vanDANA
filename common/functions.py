from dolfin import *
from scipy.interpolate import splev
import numpy as np
import sys, math, os, cppimport
from ufl import Jacobian  

# Functions for boundary conditions
class Inflow(UserExpression):

    def __init__(self, param, mesh, **kwargs):
        super(Inflow, self).__init__(**kwargs)
        self.param = param
        self.mesh = mesh

    def eval_cell(self, values, x, ufc_cell):

        # Create DOLFIN Cell
        cell = Cell(self.mesh, ufc_cell.index)
        
        # Get normal for current facet
        assert(ufc_cell.local_facet >= 0)
        n = cell.normal(ufc_cell.local_facet)
        
        # Compute boundary value
        t = self.param["time"].t
        period = self.param["period"]
        nm = self.param["nm"].cycle
        Area = self.param["Area"]
        Vsc = self.param["Vsc"]
        Tsc = self.param["Tsc"]
        func = self.param["func"]

        val = (splev(t*Tsc - nm*period*Tsc, func)/Area)/Vsc
        values[0] = -n.x()*val
        values[1] = -n.y()*val
        values[2] = -n.z()*val 

    def value_shape(self):
        return (3,)


code = '''
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include <dolfin/function/Expression.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/geometry/Point.h>

class Inflow_x : public dolfin::Expression
{
public:
    
    double v;
    std::shared_ptr<dolfin::MeshFunction<std::size_t>> cell_data;

    Inflow_x(double v_, std::shared_ptr<dolfin::MeshFunction<std::size_t>> cell_data_) : Expression() {
        v = v_;
        cell_data = cell_data_;
    }

    void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x, const ufc::cell& c) const override
    {
        assert(cell_data);
        const dolfin::Cell cell(*cell_data->mesh(), c.index);
        
        assert(c.local_facet >= 0);
        dolfin::Point n = cell.normal(c.local_facet);
        values[0] = -n.x()*v;  
    }
};

class Inflow_y : public dolfin::Expression
{
public:
    
    double v;
    std::shared_ptr<dolfin::MeshFunction<std::size_t>> cell_data;

    Inflow_y(double v_, std::shared_ptr<dolfin::MeshFunction<std::size_t>> cell_data_) : Expression() {
        v = v_;
        cell_data = cell_data_;
    }

    void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x, const ufc::cell& c) const override
    {
        assert(cell_data);
        const dolfin::Cell cell(*cell_data->mesh(), c.index);
        
        assert(c.local_facet >= 0);
        dolfin::Point n = cell.normal(c.local_facet);
        values[0] = -n.y()*v;  
    }
};

class Inflow_z : public dolfin::Expression
{
public:
    
    double v;
    std::shared_ptr<dolfin::MeshFunction<std::size_t>> cell_data;

    Inflow_z(double v_, std::shared_ptr<dolfin::MeshFunction<std::size_t>> cell_data_) : Expression() {
        v = v_;
        cell_data = cell_data_;
    }

    void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x, const ufc::cell& c) const override
    {
        assert(cell_data);
        const dolfin::Cell cell(*cell_data->mesh(), c.index);
        
        assert(c.local_facet >= 0);
        dolfin::Point n = cell.normal(c.local_facet);
        values[0] = -n.z()*v;  
    }
};

PYBIND11_MODULE(SIGNATURE, m)
{
    py::class_<Inflow_x, std::shared_ptr<Inflow_x>, dolfin::Expression>(m, "Inflow_x")
    .def(py::init<double, std::shared_ptr<dolfin::MeshFunction<std::size_t>>>())
    .def_readwrite("v", &Inflow_x::v)
    .def_readwrite("cell_data", &Inflow_x::cell_data);
    
    py::class_<Inflow_y, std::shared_ptr<Inflow_y>, dolfin::Expression>(m, "Inflow_y")
    .def(py::init<double, std::shared_ptr<dolfin::MeshFunction<std::size_t>>>())
    .def_readwrite("v", &Inflow_y::v)
    .def_readwrite("cell_data", &Inflow_y::cell_data);
    
    py::class_<Inflow_z, std::shared_ptr<Inflow_z>, dolfin::Expression>(m, "Inflow_z")
    .def(py::init<double, std::shared_ptr<dolfin::MeshFunction<std::size_t>>>())
    .def_readwrite("v", &Inflow_z::v)
    .def_readwrite("cell_data", &Inflow_z::cell_data);   
}
'''



class Outflow(UserExpression):

    def __init__(self, param, mesh, **kwargs):
        super(Outflow, self).__init__(**kwargs)
        self.param = param
        self.mesh = mesh

    def eval_cell(self, values, x, ufc_cell):

        # Create DOLFIN Cell
        cell = Cell(self.mesh, ufc_cell.index)
        
        # Get normal for current facet
        assert(ufc_cell.local_facet >= 0)
        n = cell.normal(ufc_cell.local_facet)
        
        # Compute boundary value
        t = self.param["time"].t
        period = self.param["period"]
        nm = self.param["nm"].cycle
        Area = self.param["Area"]
        Vsc = self.param["Vsc"]
        Tsc = self.param["Tsc"] 
        func = self.param["func"]

        val = splev(t*Tsc - nm*period*Tsc, func)*133.322/(1060*Vsc*Vsc)
        values[0] = val

    def value_shape(self):
        return ()

class Temperature_balloon(UserExpression):

    def __init__(self, param, mesh, **kwargs):
        super(Temperature_balloon, self).__init__(**kwargs)
        self.param = param
        self.mesh = mesh

    def eval_cell(self, values, x,  ufc_cell):

        # Create DOLFIN Cell
        cell = Cell(self.mesh, ufc_cell.index)

        # Get normal for current facet
        assert(ufc_cell.local_facet >= 0)
        n = cell.normal(ufc_cell.local_facet)

        # Compute boundary value
        t = self.param["time"].t
        Tsc = self.param["Tsc"]
        func = self.param["func"]

        val = splev(t*Tsc, func)
        values[0] = val

    def value_shape(self):
        return ()


# Algorithm for remeshing / improving mesh quality
def mesh_smoothening(mesh):

    dim = mesh.geometry().dim()
    x = SpatialCoordinate(mesh)
    dx_dxi = Jacobian(mesh)
    Jc = abs(det(dx_dxi))
        
    def grad_xi(f):
        df_dx = grad(f)
        return dot(df_dx,dx_dxi)
        
    h = mesh.hmax()
    V = VectorFunctionSpace(mesh,"CG", 1)

    # Variational problem for smoothing:
    u = TrialFunction(V)
    v = TestFunction(V)
    x_smoothed = x + u

    # Parameter to control smoothing effect;
    # larger parameter => stronger smoothing:
    smoothing_strength = Constant(1e4)

    # Penalize both large changes of position over a single element and deviation
    # from the unsmoothed position, with smoothing_strength deciding the relative
    # weights.
    dxz = Measure("dx", metadata={"quadrature_degree":2}, domain=mesh)
    res = (smoothing_strength*inner(grad_xi(x_smoothed),grad_xi(v))
           + dot(u,v))*(1/Jc)*dxz

    # Solve for displacement to deform mesh; this does require a linear solve,
    # but it should be efficient to approximate using an iterative solver
    # with a loose tolerance.
    uh = Function(V)
    solve(lhs(res)==rhs(res),uh,
          bcs=[DirichletBC(V,Constant((0,)*dim),"on_boundary")],
          solver_parameters={"linear_solver" : "cg",
                                 "krylov_solver":{"relative_tolerance":1e-4}})

    # Deform by displacement from variational problem and plot smoothed mesh:
    ALE.move(mesh, project(uh, V))
        
    ratios = MeshQuality.radius_ratio_min_max(mesh)
    ratio_min = ratios[0]
    ratio_max = ratios[1]

    return mesh, ratio_min, ratio_max



# Function to create nullspace
def attach_nullspace(Ap, x_, Q):

    """Create null space basis object and attach to Krylov solver."""
    null_vec = Vector(x_.vector())
    Q.dofmap().set(null_vec, 1.0)
    null_vec *= 1.0 / null_vec.norm('l2')
    Aa = as_backend_type(Ap)
    null_space = VectorSpaceBasis([null_vec])
    Aa.set_nullspace(null_space)
    Aa.null_space = null_space
    return null_space




# Function to calculate denominator for courant number
def DENO(u, u_components, Mpi, mesh, h_f_X):

    DN_local = 0
    NM_local = 0
    vertex_values_h_f_X = h_f_X.compute_vertex_values(mesh)
    vertex_mag_u = np.zeros(len(vertex_values_h_f_X))

    for ui in range(u_components):
        vertex_values_u = u[ui].compute_vertex_values(mesh) 
        DN_local += np.max(np.abs(vertex_values_u / vertex_values_h_f_X))
        vertex_mag_u += np.square(vertex_values_u)

    NM_local =  np.max(np.sqrt(vertex_mag_u)*vertex_values_h_f_X)

    DN = Mpi.Max(DN_local)
    NM = Mpi.Max(NM_local)

    return DN, NM


# Function to assign-vector in MPI
def vector_assign_in_parallel(v, w):

    v.vector().set_local(w.vector().get_local()[:])
    v.vector().apply("insert")


# Degrees of freedom
def Calc_total_DOF(Mpi, **kwargs):

    DOFS = dict()
    for key, value in kwargs.items():
        
        dof = 0
        for i in value:
            dof += Mpi.Sum(i.vector().get_local().size)
        DOFS[key] = dof

    return DOFS


# Delta-interpolation
def interpolate_nonmatching_mesh_delta(fsi_interpolation, u0, V, abc, flag):

    u = Function(V)
    u.vector().zero()
    u0.set_allow_extrapolation(True)
    fsi_interpolation.interpolate_delta(u0, u, abc, flag)
    return u 


# Other miscellaneous functions

# Divergence of a vector
def divergence(u, u_components):
    
    DIV = 0
    for ui in range(u_components): DIV += u[ui].dx(ui) 
    return DIV

# Viscous dissipation source term (for energy equation)
def Qf(u, Ec, Re):
    
    return inner(((2*Ec)/Re)*sym(nabla_grad(u)), nabla_grad(u))

# Perfusion equation (quadratic)
def PFE(Tf_n): 

    return conditional(Tf_n >= 0.725, (6.18*Tf_n*Tf_n) - (7.39*Tf_n) + 2.21, 0.1)  

# Returns a value rounded down to a specific number of decimal places
def round_decimals_down(number:float, decimals:int=8):

    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more")
    elif decimals == 0:
        return math.floor(number)

    factor = 10 ** decimals
    return math.floor(number * factor) / factor

# Calculate and print runtime statistics and update timestep
def calc_runtime_stats_timestep(Mpi, u, u_components, t, tsp, text_file_handles, mesh, hmin_f, h_f_X, Re, time_control): 

    DN, NM = DENO(u, u_components, Mpi, mesh, h_f_X)
    C_no_real = DN*tsp
    local_Re = NM*float(Re) 
    tsp_min = 0.125*hmin_f

    Mpi.set_barrier()
    if Mpi.get_rank() == 0:
        text_file_handles[1].write(f"{t}    {tsp}     {C_no_real}     {local_Re}\n")

    if time_control['adjustable_timestep'] == True:

        if C_no_real > time_control['C_no']:
            tsp = round_decimals_down((0.25*time_control['C_no'] + 0.75*C_no_real)/DN, 5)
        if tsp <= tsp_min:
            tsp = round_decimals_down(tsp_min, 5)     

    return tsp 

# Update solution at end of time loop
def update_variables(update, u_components, problem_physics):

    u_ = update[0]
    p_ = update[1]

    for ui in range(u_components):
        u_[2][ui].assign(u_[1][ui])
        u_[1][ui].assign(u_[0][ui]) 
    p_[2].assign(p_[1])    
    p_[1].assign(p_[0])
    
    T_ = update[2]

    T_[2].assign(T_[1])
    T_[1].assign(T_[0])

    if problem_physics['solve_FSI'] == True:

        Dp_ = update[3]
        Lm_ = update[4]

        Dp_[2].assign(Dp_[1])
        Lm_[1].assign(Lm_[0])

        if problem_physics['solve_temperature'] == True:

            Ts_   = update[5]
            LmTs_ = update[6]    
             
            Ts_[1].assign(Ts_[0])
            LmTs_[1].assign(LmTs_[0])
            



