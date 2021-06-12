import numpy as np


def sphe2cart(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z


def cart2sphe(x, y, z):
    r = np.linalg.norm((x, y, z))
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    return r, theta, phi


def isect_line_plane_v3_4d(p0, p1, plane, epsilon=1e-6):
    """
    p0, p1: define the line
    p_co, p_no: define the plane:
        p_co is a point on the plane (plane coordinate).
        p_no is a normal vector defining the plane direction;
             (does not need to be normalized).

    return a Vector or None (when the intersection can't be found).
    """
    p_no = plane[:3]
    u = p1 - p0
    dot = np.dot(p_no, u)
    if dot > epsilon:
        # calculate a point on the plane
        # (divide can be omitted for unit hessian-normal form).
        p_co = p_no * (-plane[3] / np.linalg.norm(p_no))

        w = p0 - p_co
        fac = -np.dot(p_no, w) / dot
        u = u * fac
        return p0 + u
    else:
        return None


def find_plane_ray_intersection(hull, centroid, ray):
    dist = float('inf')
    intersection = []
    for plane in hull.equations:
        inte = isect_line_plane_v3_4d(centroid, ray, plane)
        if inte is not None:
            new_dist = np.linalg.norm(inte - centroid)
            if new_dist < dist:
                dist = new_dist
                intersection = inte
    return np.array(intersection)


def get_coord_point_list(hull, centroid, point, total_layers,
                         initial_distance, delta_r):
    intersect = centroid - find_plane_ray_intersection(hull, centroid, point)
    r, theta, phi = cart2sphe(*intersect)
    r += initial_distance
    all_points = []
    for layer in range(total_layers):
        all_points.append(centroid - (sphe2cart(r, theta, phi)))
        r += delta_r
    return all_points


def generate_points(hull, step, initial_distance, total_layers, delta_r):
    # step = 30
    # initial_distance = 2.5
    # total_layers
    # delta_r
    centroid = np.mean(hull.points[hull.vertices, :], axis=0)

    cur_point = centroid - sphe2cart(1, 0, 0)
    new_cur_point = centroid - cur_point
    r, theta, phi = cart2sphe(*new_cur_point)

    all_points = get_coord_point_list(hull, centroid, cur_point, total_layers,
                                      initial_distance, delta_r)

    for theta in range(step, 180, step):
        for phi in range(0, 360, step):
            # move center of spherical coordinate to center of convex hull
            cur_point = centroid - sphe2cart(1,
                                             np.radians(theta),
                                             np.radians(phi))
            all_points += get_coord_point_list(hull, centroid, cur_point,
                                               total_layers, initial_distance,
                                               delta_r)

    return np.array(all_points)
