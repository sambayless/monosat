################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_UPPER_SRCS += \
../minilib/Minisat_Solver.C 

OBJS += \
./minilib/Minisat_Solver.o 

C_UPPER_DEPS += \
./minilib/Minisat_Solver.d 


# Each subdirectory must supply rules for building sources it contributes
minilib/%.o: ../minilib/%.C
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -I"/home/sam/workspaceC/modsat" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


