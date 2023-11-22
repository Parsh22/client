import re
from functools import reduce

amino_acids_dict = {}
order_evals_dict = {}
decoded_words_dict = {}
codon_sequences = ()


def expand_sequence(accumulator, character):
    if character == "{":
        return accumulator
    elif character == "}":
        return accumulator
    elif character.isdigit():
        return accumulator[:-1] + accumulator[-1] * (int(character) - 1)
    else:
        return accumulator + character


def read_codons(codon_file):
    amino_acids_dict.clear()

    file_obj = open(codon_file)

    amino_acid_regex = re.compile(r"([A-Z][a-zA-Z]*): (.+)")

    for line in file_obj:
        invalid_sequence = False

        match = amino_acid_regex.match(line)

        if not match:
            continue

        if not invalid_sequence:
            amino_acid = match.group(1)
            sequences = match.group(2)

            sequence_array = []

            for seq in sequences.split(", "):
                if not invalid_sequence:
                    expanded_sequence = ""
                    index = 0
                    while index < len(seq):
                        if seq[index] == "{":
                            number = ""
                            index += 1
                            while seq[index] != "}":
                                number += seq[index]
                                index += 1
                            expanded_sequence += expanded_sequence[-1] * (
                                int(number) - 1
                            )
                            index += 1
                        else:
                            expanded_sequence += seq[index]
                            index += 1

                    expanded_sequence = expanded_sequence.rstrip()

                    if re.search(r"[^AUGC]", expanded_sequence):
                        invalid_sequence = True
                        continue

                    sequence_array.append(expanded_sequence)

                    if sequence_array:
                        amino_acids_dict[amino_acid] = sequence_array

    return amino_acids_dict


def read_evals(eval_file):
    order_evals_dict.clear()

    file_obj = open(eval_file)
    order_regex = re.compile(r"([A-Za-z0-9]+): (L|R), (PO|PR|I)")

    for line in file_obj:
        evaluation = re.fullmatch(order_regex, line.strip())
        if evaluation:
            order_evals_dict[evaluation.group(1)] = [
                evaluation.group(2),
                evaluation.group(3),
            ]

    return None


def find_max_length_codon(accumulator, current):
    return accumulator if len(accumulator) >= len(current) else current


def encode(sequence):
    encoded_sequence = ""
    for codon in sequence.split(" "):
        if codon in amino_acids_dict:
            seq = amino_acids_dict[codon]

            if len(seq) > 1:
                max_length_codon = reduce(find_max_length_codon, seq)
                encoded_sequence += max_length_codon

            else:
                encoded_sequence += seq[0]
        else:
            continue

    return encoded_sequence


def decode(sequence, associations=None):
    index = 0
    decoded_sequence = ""
    while index < len(sequence):
        max_length = 0
        amino_acid_match = ""
        substring = ""

        for amino_acid, sequences in amino_acids_dict.items():
            for seq in sequences:
                if sequence[index : index + len(seq)] == seq:
                    if len(seq) > max_length:
                        max_length = len(seq)
                        amino_acid_match = amino_acid
                        substring = sequence[index : index + len(seq)]
        if max_length > 0:
            decoded_sequence += amino_acid_match + " "
            if associations != None:
                associations.append((amino_acid_match, substring))
            index += max_length
        else:
            index += 1
    return decoded_sequence.strip()


def helper_exchange(tuple_item):
    array = amino_acids_dict[tuple_item[0]]

    for item in array:
        if item != tuple_item[1]:
            return (tuple_item[0], item)

    return tuple_item


def operate(sequence, eval_name):
    if eval_name not in order_evals_dict:
        return None

    eval_information = order_evals_dict[eval_name]

    print("Original sequence is: ", sequence)

    final_operations = []

    operations = ["START", "STOP", "SWAP", "EXCHANGE", "DEL"]

    operations_count = 0

    direction = eval_information[0]

    order = eval_information[1]

    codon_sequences = []

    if direction == "R":
        sequence = sequence[::-1]

    decoded = decode(sequence, codon_sequences)

    print("Decoded is: ", decoded)

    print("Sequences is right now: ", codon_sequences)

    adding_section = False

    for pair in codon_sequences:
        if pair[0] == "START":
            adding_section = True
        elif pair[0] == "STOP":
            adding_section = False
        elif adding_section:
            final_operations.append(pair)

    final_operations.insert(0, ("START", amino_acids_dict["START"][0]))

    final_operations.append(("STOP", amino_acids_dict["STOP"][0]))

    codon_sequences = final_operations

    adding_section = False
    final_operations = []
    index = 0

    while index < len(codon_sequences):
        pair = codon_sequences[index]
        word = codon_sequences[index][0]

        if order == "PO":
            # Postfix logic
            if word == "START":
                adding_section = True
                final_operations.append(pair)

            elif word == "STOP":
                adding_section = False
                final_operations.append(pair)

            elif word != "START" and word != "STOP" and adding_section:
                if word == "SWAP":
                    if len(final_operations) >= 3:
                        if (
                            final_operations[len(final_operations) - 1][0]
                            not in operations
                            and final_operations[len(final_operations) - 2][0]
                            not in operations
                        ):
                            # Swap the last two elements
                            last_two = final_operations[len(final_operations) - 1]
                            final_operations[
                                len(final_operations) - 1
                            ] = final_operations[len(final_operations) - 2]
                            final_operations[len(final_operations) - 2] = last_two

                            # Increment operation count
                            operations_count += 1

                        elif (
                            final_operations[len(final_operations) - 1][0] != "START"
                            and final_operations[len(final_operations) - 2][0]
                            != "START"
                        ):
                            final_operations.append(pair)
                            operations_count += 1

                elif word == "DEL":
                    # Remove the last element
                    if len(final_operations) >= 2 and index > 1:
                        if (
                            final_operations[len(final_operations) - 1][0]
                            not in operations
                        ):
                            final_operations.pop(len(final_operations) - 1)
                            operations_count += 1
                        elif final_operations[len(final_operations) - 1][0] != "START":
                            final_operations.append(pair)
                            operations_count += 1

                elif word == "EXCHANGE":
                    if len(final_operations) >= 2 and index > 1:
                        if (
                            final_operations[len(final_operations) - 1][0]
                            not in operations
                        ):
                            final_operations[
                                len(final_operations) - 1
                            ] = helper_exchange(
                                final_operations[len(final_operations) - 1]
                            )
                            operations_count += 1

                        elif final_operations[len(final_operations) - 1][0] != "START":
                            final_operations.append(pair)
                            operations_count += 1

                else:
                    final_operations.append(pair)

        elif order == "PR":
            # Prefix logic
            if word == "START":
                print("Current tuple:", final_operations)
                adding_section = True
                final_operations.append(pair)

            elif word == "STOP":
                print("Current tuple:", final_operations)
                adding_section = False
                final_operations.append(pair)

            elif adding_section and word != "START" and word != "STOP":
                if word == "SWAP":
                    print("Current tuple:", final_operations)
                    if index < len(codon_sequences) - 2:
                        if (
                            codon_sequences[index + 1][0] not in operations
                            and codon_sequences[index + 2][0] not in operations
                        ):
                            print("In swap for PR")
                            final_operations.append(codon_sequences[index + 2])
                            final_operations.append(codon_sequences[index + 1])
                            index += 2

                            operations_count += 1

                        elif (
                            codon_sequences[index + 1][0] != "STOP"
                            and codon_sequences[index + 2][0] != "STOP"
                        ):
                            final_operations.append(pair)
                            operations_count += 1

                elif word == "DEL":
                    print("Current tuple:", final_operations)
                    if index < len(codon_sequences) - 1:
                        if codon_sequences[index + 1][0] not in operations:
                            index += 1
                            operations_count += 1
                        elif codon_sequences[index + 1][0] != "STOP":
                            final_operations.append(pair)
                            operations_count += 1

                elif word == "EXCHANGE":
                    if index < len(codon_sequences) - 1:
                        if codon_sequences[index + 1][0] not in operations:
                            print(
                                "Sequences at index + 1: ", codon_sequences[index + 1]
                            )
                            print(
                                "Final operations before exchange: ", final_operations
                            )
                            final_operations.append(
                                helper_exchange(codon_sequences[index + 1])
                            )
                            print("Final operations after exchange: ", final_operations)
                            index += 1
                            operations_count += 1

                        elif codon_sequences[index + 1][0] != "STOP":
                            final_operations.append(pair)
                            operations_count += 1

                else:
                    print("Current tuple:", final_operations)
                    decoded_array = amino_acids_dict[word]

                    final_operations.append(pair)

        elif order == "I":
            # Infix logic
            if word == "START":
                adding_section = True
                final_operations.append(pair)

            elif word == "STOP":
                adding_section = False
                final_operations.append(pair)

            elif word != "START" and word != "STOP" and adding_section:
                if word == "SWAP":
                    # Swap only if its the third one from beginning and end
                    print("Index is: ", index)
                    if index < len(codon_sequences) - 1 and index > 1:
                        if (
                            final_operations[len(final_operations) - 1][0]
                            not in operations
                            and codon_sequences[index + 1][0] not in operations
                        ):
                            final_operations.insert(
                                len(final_operations) - 1, codon_sequences[index + 1]
                            )
                            index += 1
                            operations_count += 1

                        elif (
                            final_operations[len(final_operations) - 1][0]
                            not in operations
                            and codon_sequences[index + 1][0] != "STOP"
                        ):
                            final_operations.append(pair)
                            operations_count += 1

                elif word == "DEL":
                    if index < len(codon_sequences) - 1:
                        if codon_sequences[index + 1][0] not in operations:
                            index += 1
                            operations_count += 1
                        elif codon_sequences[index + 1][0] != "STOP":
                            final_operations.append(pair)
                            operations_count += 1

                elif word == "EXCHANGE":
                    if index < len(codon_sequences) - 1:
                        if codon_sequences[index + 1][0] not in operations:
                            final_operations.append(
                                helper_exchange(codon_sequences[index + 1])
                            )
                            index += 1
                            operations_count += 1

                        elif codon_sequences[index + 1][0] != "STOP":
                            final_operations.append(pair)
                            operations_count += 1

                else:
                    final_operations.append(pair)

        index += 1

    if operations_count == 0:
        simplified = []

        for pair in final_operations:
            if pair[0] not in operations:
                simplified.append(pair)

        final_operations = simplified

        result = ""
        for element in final_operations:
            result += element[1]

        final_operations = result

    elif operations_count > 0:
        concatenated = ""

        for element in final_operations:
            concatenated += element[1]

        if direction == "R":
            concatenated = concatenated[::-1]

        final_operations = operate(concatenated, eval_name)

    return final_operations
