clear
for i in {100,500,1000,10000}
do
 echo "$i"
done
k=100

gcc 2DPoiseuille2.cpp -o P
./P "$k"

echo "done"