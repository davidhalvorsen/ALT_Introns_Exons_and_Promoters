#!/usr/bin/env python
import random
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

    current_list_of_telomere_lengths = [telomere-50 for telomere in current_list_of_telomere_lengths]
    current_division_number += 1
