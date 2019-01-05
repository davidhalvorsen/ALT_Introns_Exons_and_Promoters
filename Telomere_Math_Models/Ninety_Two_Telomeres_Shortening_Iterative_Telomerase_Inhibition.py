#!/usr/bin/env python
import random
# 6*767 is more than required for immortality
# telomerase_bp_added_per_division = 6*767
# 383 is less than half of what's required for immortality
# telomerase_bp_added_per_division = 6*383
telomerase_bp_added_per_division = 6*383
# telomerase_bp_added_per_division = 0
# telomerase_stopped_at_division_number = 600
telomerase_stopped_at_division_number = 5
shortest = 5000
longest = 15000
num_of_telomeres = 92
telomeres_longer_than_3000 = True
current_division_number = 1

def initialze_92_telomere_lengths(shortest, longest, num_of_telomeres):
    telomere_length_list = []
    for j in range(num_of_telomeres):
        telomere_length_list.append(random.randint(shortest, longest))
    return telomere_length_list
#print "This is the initial list of telomere lengths: "
initial_list_of_telomere_lengths = initialze_92_telomere_lengths(shortest, longest, num_of_telomeres)
#print(initial_list_of_telomere_lengths)
current_list_of_telomere_lengths = initial_list_of_telomere_lengths

while telomeres_longer_than_3000 == True:
    # current division
    print "It is division # " + str(current_division_number)
    print "This is the current list of telomere lengths: "
    print current_list_of_telomere_lengths

    current_telomere_number = 1
    for telomere in current_list_of_telomere_lengths:
        if telomere <= 5000 and telomere > 3000:
            print "Telomere number " + str(current_telomere_number) + " is only " + str(telomere) + " bp long."
            print "p53-mediated senescence would be triggered, BUT p53 is mutated!"
        elif telomere <= 3000:
            print "Telomere number " + str(current_telomere_number) + " is only " + str(telomere) + " bp long."
            print "There is massive genomic instability and cell death!"
            telomeres_longer_than_3000 = False
            break
        current_telomere_number += 1
    # adding bp with telomerase
    telomere_bp_left_to_add = telomerase_bp_added_per_division
    if telomerase_stopped_at_division_number <= current_division_number:
        telomere_bp_left_to_add = 0
    while telomere_bp_left_to_add > 0:
        # get index for shortest telomere
        shortest_telomere_index = current_list_of_telomere_lengths.index(min(current_list_of_telomere_lengths))
        # get length of shortest telomere
        length_of_shortest_telomere = current_list_of_telomere_lengths[shortest_telomere_index]
        # add 6bp to shortest telomere
        new_length_of_shortest_telomere = current_list_of_telomere_lengths[shortest_telomere_index] + 6
        # change list entry at index # to new telomere length
        current_list_of_telomere_lengths[shortest_telomere_index]=new_length_of_shortest_telomere
        # less telomerase left
        telomere_bp_left_to_add = telomere_bp_left_to_add - 6
    # setting up for next division
    current_list_of_telomere_lengths = [telomere-50 for telomere in current_list_of_telomere_lengths]
    current_division_number += 1
