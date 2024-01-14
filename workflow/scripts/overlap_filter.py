# Usage: python3 overlap_filter.py <input>.bed <output>.bed
from sys import argv


input_bed = argv[1]

# init indexes
new_bed_list = []
prev_start = 0
prev_end = 0
prev_score = 0
prev_length = 0

with open(input_bed) as bed:
    for line in bed:

        if line[:5] == 'track':
            continue

    	# read line
        hor = line.strip().split('\t')
        start = int(hor[1])
        end = int(hor[2])
        length = end - start
        score = float(hor[4])
        # check if starts or ends match
        if start in range(prev_start-10, prev_start+10) or end in range(prev_end-10, prev_end+10):
            if prev_score-0.01 < score < prev_score+0.01:
                #print(start, prev_score, score)
                # scores are same, look at length
                if length > prev_length:
                    new_bed_list.pop()
                else:
                    continue
            elif score > prev_score:
            	# delete previous hor if list is not empy
            	if len(new_bed_list) > 0:
                    new_bed_list.pop()
            else:
            	# don't add current hor
                continue
        # add hor
        new_bed_list.append(line)
        # update indexes
        prev_start = start
        prev_end = end
        prev_score = score
        prev_length = length

for hor in new_bed_list:
    print(hor.strip())
