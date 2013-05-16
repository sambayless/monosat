################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../simp/Main.o \
../simp/SimpSolver.o 

CC_SRCS += \
../simp/Main.cc \
../simp/SimpSolver.cc 

OBJS += \
./simp/Main.o \
./simp/SimpSolver.o 

CC_DEPS += \
./simp/Main.d \
./simp/SimpSolver.d 


# Each subdirectory must supply rules for building sources it contributes
simp/%.o: ../simp/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -DINTERPOLATION -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


