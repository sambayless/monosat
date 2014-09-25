################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../aiger/aiger.c 

OBJS += \
./aiger/aiger.o 

C_DEPS += \
./aiger/aiger.d 


# Each subdirectory must supply rules for building sources it contributes
aiger/%.o: ../aiger/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -DNDEBUG -DCGAL_DISABLE_ROUNDING_MATH_CHECK -I.././ -I.././ -O3 -Wall -c -fmessage-length=0 -frounding-math -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


