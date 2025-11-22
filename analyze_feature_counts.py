with open("output/Normal_2_featureCounts.txt", "r") as file:
    for idx,line in enumerate(file):
        # print(line.split("\t")[0])
        print(line)
        print("----")
        if idx == 5:
            break