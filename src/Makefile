OS=win
# OS=unix

F90=gfortran
# F90=ifort

DO=./
DB=./
DS=./
F=vxc2data
EXE=$(F)

ifeq ($(F90),gfortran) 

    Csw = -c  -O2 -finit-local-zero  -fbounds-check -ffree-line-length-none 
    Lsw = 

    OBJ = $(DO)write_vtk.o \
    $(DO)utilites.o \
    $(DO)$(F).o 

    EXE: $(OBJ)
		$(F90) $(OBJ) $(Lsw)  -o $(DB)$(EXE) 

    $(OBJ): $(DO)%.o: $(DS)%.f90
		$(F90) $(Csw) $< -o $@
endif

ifeq ($(OS),unix) 
    ifeq ($(F90),ifort) 
        Csw = -c -O2 -zero 
    Lsw = 

    OBJ = $(DO)write_vtk.obj \
    $(DO)utilites.obj \
    $(DO)$(F).obj

    EXE: $(OBJ)
		$(F90) $(OBJ) $(Lsw)  -o $(DB)$(EXE) 

    $(OBJ): $(DO)%.obj: $(DS)%.f90
		$(F90) $(Csw) $< -o $@

    endif
endif

ifeq ($(OS),win) 
    ifeq ($(F90),ifort) 
        Csw = -c /QO2  /check:bounds /Qzero  /heap-arrays:2000 /assume:buffered_io
        # /heap-arrays:2000
        Lsw = 

        OBJ = $(DO)write_vtk.obj \
        $(DO)utilites.obj \
        $(DO)$(F).obj

        EXE: $(OBJ)
			$(F90) $(OBJ) $(Lsw)  -o $(DB)$(EXE) 

        $(OBJ): $(DO)%.obj: $(DS)%.f90
			$(F90) $(Csw) $< -o $@

    endif
endif


# Cleaning everything
clean:
ifeq ($(OS),win) 
    ifeq ($(F90),gfortran) 
		del *.o
    endif
    ifeq ($(F90),ifort) 
		del *.obj
    endif
cleanvtk:
	del *.vtk
endif

ifeq ($(OS),unix) 
    ifeq ($(F90),gfortran) 
		rm *.o
    endif

    ifeq ($(F90),ifort) 
		rm *.obj
    endif

cleanvtk:
	rm *.vtk
endif
