scores=[100,80,40,5,1]
best_k=-1
best_total=-1
for k in range(32):
    count=0
    total=0
    for i in range(5):
        bit=((k>>i)&1)
        count=count+bit
        total=total+bit*scores[i]
    if(count==3):
        if(total>best_total):
            best_total=total
            best_k=k
print("The best total we get is ",best_total,", and the best k we have such that it accomodates 3 is",best_k)