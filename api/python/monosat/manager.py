from monosat.monosat_c import Monosat
#Like singleton, but one instance per _solver_, not per process
class Manager(type):

    def __call__(cls, *args, **kwargs):        
        if cls not in Monosat()._getManagers():
            inst = super(Manager, cls).__call__(*args, **kwargs)
            inst._solver=Monosat().getSolver()
            Monosat()._getManagers()[cls] = inst            
        return Monosat()._getManagers()[cls]