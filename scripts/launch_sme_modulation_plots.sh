
observable="n_bjets"

for year in 2016 2017
do
    for sample in signal singletop
    do
        for wilson in cLXX cLXY cLXZ cLYZ cRXX cRXY cRXZ cRYZ cXX cXY cXZ cYZ dXX dXY dXZ dYZ
        do
            python bin/sme_time.py ${observable} ${year} ${wilson} ${sample}
        done
    done
done
