# Initialize all three methods
# for method in ST DOP MCS
for method in MCS
do
    python init3.py -l 32 1 3 $method
    python init3.py -l 64 1 4 $method
    python init3.py -l 64 2 4 $method
done

# 32-run res III, m= 1
for i in {4..20}
do 
    python extend3_st.py -l 32 1 $i 3
    python extend3_dop.py -l 32 1 $i 3
done

# 32-run res III, m= 1
for i in {2..20}
do 
    python extend3_mcs.py -l 32 1 $i 3
done


# 64-run, res IV, m = 1
for i in {5..15} 
do 
    python extend3_st.py -l 64 1 $i 4
    python extend3_dop.py -l 64 1 $i 4
done

# 64-run, res IV, m = 1
for i in {3..15} 
do 
    python extend3_mcs.py -l 64 1 $i 4
done

# 64-run, res IV, m = 2
for i in {3..12} 
do 
    python extend3_st.py -l 64 2 $i 4
    python extend3_dop.py -l 64 2 $i 4
done

# 64-run, res IV, m = 2
for i in {2..12} 
do 
    python extend3_mcs.py -l 64 2 $i 4
done
