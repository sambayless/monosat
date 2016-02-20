################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../dgl/alg/DisjointSets.cpp \
../dgl/alg/LinkCutCost.cpp 

OBJS += \
./dgl/alg/DisjointSets.o \
./dgl/alg/LinkCutCost.o 

CPP_DEPS += \
./dgl/alg/DisjointSets.d \
./dgl/alg/LinkCutCost.d 


# Each subdirectory must supply rules for building sources it contributes
dgl/alg/%.o: ../dgl/alg/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS -DNDEBUG -I.././ -O3 -Wall -c -fmessage-length=0 -std=c++11 -Wno-unused-variable -Wno-unused-but-set-variable  -static  -static-libgcc     -static-libstdc++ -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


