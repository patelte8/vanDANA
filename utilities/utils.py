from dolfin import MPI, Timer, Constant
import sys, subprocess, os
from os import getpid


# Printing helper functions 
RED = "\033[1;37;31m%s\033[0m"
BLUE = "\033[1;37;34m%s\033[0m"
GREEN = "\033[1;37;32m%s\033[0m"



# Disable printing
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore printing
def enablePrint():
    sys.stdout = sys.__stdout__	



# MPI utility
class MPI_Manage:

    def __init__(self):
        
        self.mpi_comm = MPI.comm_world
        self.my_rank = MPI.rank(self.mpi_comm)

    def get_rank(self):
        
        return self.my_rank

    def set_barrier(self):
        
        MPI.barrier(self.mpi_comm)

    def get_communicator(self): 

        return self.mpi_comm

    def Max(self, x):

        return MPI.max(self.mpi_comm, x)

    def Min(self, x):

        return MPI.min(self.mpi_comm, x)        

    def Sum(self, x):
        
        return MPI.sum(self.mpi_comm, x)       




# Counter functions in loop
def create_counters(i):

    my_list = []
    for x in range(0, i):
        x = 0; my_list.append(x)

    return my_list   

def reset_counter(j, *args):

    if args:
        for x in args:
            j[x] = 0
    else:
        for i in range(len(j)):
            j[i] = 0    

def update_counter(j, *args):
    if args:
         for x in args:
            j[x] += 1         
    else:
        for i in range(len(j)):
            j[i] += 1    



# Memory functions
def getMemoryUsage(rss=True):
    mypid = str(getpid())
    rss = "rss" if rss else "vsz"
    process = subprocess.Popen(['ps', '-o', rss, mypid],
                                stdout=subprocess.PIPE)
    out, _ = process.communicate()
    mymemory = out.split()[1]
    return eval(mymemory) / 1024

class MemoryUsage:
    def __init__(self, s):
        self.memory = 0
        self.memory_vm = 0
        self(s)

    def __call__(self, s, verbose=False):
        self.prev = self.memory
        self.prev_vm = self.memory_vm
        self.memory = MPI.sum(MPI.comm_world, getMemoryUsage())
        self.memory_vm = MPI.sum(MPI.comm_world, getMemoryUsage(False))
        if verbose:
            print(BLUE % '{0:26s}  {1:10d} MB {2:10d} MB {3:10d} MB {4:10d} MB'.format(s,
                        int(self.memory - self.prev), int(self.memory),
                        int(self.memory_vm - self.prev_vm), int(self.memory_vm)))



# Timer variables
s1 = 0.0; s2 = 0.0; s3 = 0.0; s4 = 0.0 
s5 = 0.0; s6 = 0.0; s7 = 0.0; s8 = 0.0
si = 0.0; sm = 0.0; sr = 0.0; s_dt = 0.0

timer_total = Timer("Total_run_time")
timer_dt    = Timer("Time_step_timer")
timer_s1    = Timer("Predict tentative velocity step")
timer_s2    = Timer("Pressure correction step")
timer_s3    = Timer("Velocity correction step")
timer_s4    = Timer("Energy conservation step")