

for i in $(otool -L ./pyx | tail -n 9 | head -8 | cut -f 1 -d' ' | cut -f 2 -d'/'); do
    ls $i | grep No;
done


for lib in $(ls lib*.dylib); do
    echo $lib \n ==============================
    for i in $(otool -L $lib | tail -n 9 | head -8 | cut -f 1 -d' ' | cut -f 2 -d'/'); do
        echo "#${i}-";
        ls $i | grep No;
    done
done