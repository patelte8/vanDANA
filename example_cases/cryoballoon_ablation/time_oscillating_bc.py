from dolfin import *
from scipy.interpolate import splev


# Functions to apply custom time-oscillating boundary conditions

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


xyz = '''
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include <dolfin/function/Expression.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/geometry/Point.h>

class Inflow : public dolfin::Expression
{
public:
    
    double v;
    std::shared_ptr<dolfin::MeshFunction<std::size_t>> cell_data;

    Inflow(double v_, std::shared_ptr<dolfin::MeshFunction<std::size_t>> cell_data_) : Expression(3) {
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
        values[1] = -n.y()*v;
        values[2] = -n.z()*v;
    }
};

PYBIND11_MODULE(SIGNATURE, m)
{
  py::class_<Inflow, std::shared_ptr<Inflow>, dolfin::Expression>(m, "Inflow")
    .def(py::init<double, std::shared_ptr<dolfin::MeshFunction<std::size_t>>>())
    .def_readwrite("v", &Inflow::v)
    .def_readwrite("cell_data", &Inflow::cell_data);
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