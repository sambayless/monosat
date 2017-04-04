 #/bin/sh
 
 N=20
 next_seed=1
 for M in {3,5,7}; do  
   for C in {0,10,20,30}; do
      mkdir -p  "M_${M}_C_${C}"
      for _ in {1..10}; do
	f="M_${M}_C_${C}/instance_N_${N}_M_${M}_C_${C}_${next_seed}.pcrt"
	echo $M $C $next_seed $f;      
	python3 generate_pcrt.py $N $M $C $next_seed $f
	next_seed=$((next_seed+1))
      done
   done
 done