#!/bin/bash
for i in {1..5};
do 
	python "/home/pcurran/phd/phd-scripts/run_hotspots/run_from_grid.py" "/home/pcurran/fragment-hotspot-results/patel_set/$i/out"
	
done
