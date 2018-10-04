#!/bin/bash
for i in {2..5};
do 
	python "/home/pcurran/phd/phd-scripts/run_hotspots/run_from_protein.py" "/home/pcurran/fragment-hotspot-results/patel_set/$i/protein.pdb" "-g" "/home/pcurran/src" 
	
done
