import sys
def GetDistance(point, centers):
    least_dist = 10**8   # Arbitrary large number
    for center in centers:
        dist = sum([(center[i] - point[i])**2 for i in range(len(point))])
        if dist<least_dist: 
            least_dist = dist

    return least_dist

def FarthestFirstTravel(points, k):
    centers = list()
    centers.append(points.pop(0))

    p = list(range(len(points)))    

    while len(centers)<k:
        dist_to_center = list()

        # Get the distance of the points from the centers
        for point in points: 
            dist_to_center.append(GetDistance(point, centers))

        # Get argmax of distance to centers. Add it to list of centers. Remove from points.
        argmax = max(p, key=lambda x: dist_to_center[x])
        centers.append(points.pop(argmax))

        p.pop()

    return centers

if __name__ == "__main__":
    points = list()
    flag = 1
    for point in sys.stdin:
        if point == '\n': break
        if flag:
            flag = 0
            k, m = point.split(' ')
            k, m = int(k), int(m)
            continue
        
        pt = list()
        for i in point.split(' '):
            if '\n' in i: pt.append(float(i[:-1]))
            else: pt.append(float(i))
        points.append(pt)
    
    centers = FarthestFirstTravel(points, k)

    for center in centers:
        print(' '.join([str(i) for i in center]))