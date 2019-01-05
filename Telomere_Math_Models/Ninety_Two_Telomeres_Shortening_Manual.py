#!/usr/bin/env python

import random
shortest = 5000
longest = 15000
num_of_telomeres = 92
def initialze_92_telomere_lengths(shortest, longest, num_of_telomeres):
    telomere_length_list = []
    for j in range(num_of_telomeres):
        telomere_length_list.append(random.randint(shortest, longest))
    return telomere_length_list
print "This is the initial list of telomere lengths: "
initial_list_of_telomere_lengths = initialze_92_telomere_lengths(shortest, longest, num_of_telomeres)
print(initial_list_of_telomere_lengths)

print "This is the list of telomere lengths after one division: "
list_after_one_division = [item-50 for item in initial_list_of_telomere_lengths]
print(list_after_one_division)

print "This is the list of telomere lengths after 300 divisions: "
list_after_300_divisions = [item-50*300 for item in initial_list_of_telomere_lengths]
print(list_after_300_divisions)
