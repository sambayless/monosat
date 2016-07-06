################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_UPPER_SRCS += \
../Debug.C \
../Hardware_adders.C \
../Hardware_clausify.C \
../Hardware_sorters.C \
../Main.C \
../MiniSat.C \
../PbParser.C \
../PbSolver.C \
../PbSolver_convert.C \
../PbSolver_convertAdd.C \
../PbSolver_convertBdd.C \
../PbSolver_convertSort.C \
../SatELite.C 

OBJS += \
./Debug.o \
./Hardware_adders.o \
./Hardware_clausify.o \
./Hardware_sorters.o \
./Main.o \
./MiniSat.o \
./PbParser.o \
./PbSolver.o \
./PbSolver_convert.o \
./PbSolver_convertAdd.o \
./PbSolver_convertBdd.o \
./PbSolver_convertSort.o \
./SatELite.o 

C_UPPER_DEPS += \
./Debug.d \
./Hardware_adders.d \
./Hardware_clausify.d \
./Hardware_sorters.d \
./Main.d \
./MiniSat.d \
./PbParser.d \
./PbSolver.d \
./PbSolver_convert.d \
./PbSolver_convertAdd.d \
./PbSolver_convertBdd.d \
./PbSolver_convertSort.d \
./SatELite.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.C
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -D_FILE_OFFSET_BITS=64 -I.././ADTs -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


