clear
name1=Datos
name2=Figuras

mkdir "$name1" "$name2" 
gcc 2DPoiseuille2.cpp -o P
for k in {100,500}
do 
 echo "se está calculando el paso $k"
 ./P "$k"
 echo "$k" | python CambioNombre.py
 mv *.dat "${name1}"
 rm x y u v rho
 echo "$k" | python graficas.py
 mv *.png "${name2}"
done


