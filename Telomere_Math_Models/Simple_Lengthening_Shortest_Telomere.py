#!/usr/bin/env python

# this is some testing I did to figure out how to have the while telomerase addition thing 
telomerase_bp_added_per_division = 6*767

telomere_lengths = [5583, 12162, 6389, 9970, 10351, 7970, 10257, 10570, 11226, 10003, 12418, 5145, 9841, 8593, 3683, 9251, 7884, 5070, 9813, 5316, 6170, 5067, 6420, 10141, 3465, 6754, 11450, 5883, 11551, 7127, 7170, 11200, 7788, 7079, 5492, 11992, 8606, 11573, 4600, 6108, 9872, 9659, 7811, 8597, 6907, 3161, 5425, 7053, 5638, 5540, 12154, 7962, 11386, 5722, 6675, 11306, 2967, 11208, 10478, 8560, 12603, 7522, 10607, 4044, 6816, 4172, 5127, 8953, 7438, 9860, 11952, 4076, 11505, 8247, 9124, 6027, 12625, 3031, 8991, 4777, 5177, 9458, 5578, 7392, 6019, 4637, 4912, 4458, 11358, 5594, 12343, 3874]

print telomere_lengths
# get index for shortest telomere
print telomere_lengths.index(min(telomere_lengths))

# get length of shortest telomere
print telomere_lengths[56]

# add 6bp to shortest telomere
new_length = telomere_lengths[56] + 6
print new_length

# change list entry at index # to new telomere length
telomere_lengths[56]=new_length
print telomere_lengths[56]
