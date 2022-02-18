#include "geometrycentral/surface/exact_geodesic_helpers.h"
#include "geometrycentral/surface/heat_method_distance.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vector_heat_method.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geom;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;
polyscope::CurveNetwork* psCurve;

std::vector<SurfacePoint> curve;

// helper stolen from geometrycentral/trace_geodesics.h
std::array<Vector2, 3>
vertexCoordinatesInTriangle(IntrinsicGeometryInterface& geom, Face face) {
    return {Vector2{0., 0.}, geom.halfedgeVectorsInFace[face.halfedge()],
            -geom.halfedgeVectorsInFace[face.halfedge().next().next()]};
}

// helper stolen from geometrycentral/trace_geodesics.h
Vector2 baryCoordsToFaceCoords(const std::array<Vector2, 3>& vertCoords,
                               Vector3 baryCoord) {
    return vertCoords[0] * baryCoord.x + vertCoords[1] * baryCoord.y +
           vertCoords[2] * baryCoord.z;
}

polyscope::CurveNetwork* drawCurve(ManifoldSurfaceMesh& mesh,
                                   VertexPositionGeometry& geom,
                                   const std::vector<SurfacePoint>& curve) {
    std::vector<Vector3> curvePositions;
    for (const SurfacePoint& p : curve) {
        curvePositions.push_back(p.interpolate(geom.vertexPositions));
    }
    return polyscope::registerCurveNetworkLine("source curve", curvePositions);
}

std::vector<SurfacePoint> generateCurve(ManifoldSurfaceMesh& mesh,
                                        VertexPositionGeometry& geom) {
    srand(0);
    size_t N = 50;
    // size_t iF = rand() % mesh->nFaces();
    // size_t iF = 916;
    size_t iF = 25182;
    curve.push_back(SurfacePoint(mesh.face(iF), Vector3{1, 1, 1} / 3.));
    Halfedge heCurr = curve[0].face.halfedge();
    curve.push_back(SurfacePoint(heCurr, 0.5));

    while (curve.size() + 1 < N) {
        heCurr = heCurr.twin().next();
        // if (rand() % 10 < 5) {
        if (curve.size() % 2 == 0) {
            heCurr = heCurr.next();
        }
        curve.push_back(SurfacePoint(heCurr, 0.5));
    }
    curve.push_back(SurfacePoint(heCurr.twin().face(), Vector3{1, 1, 1} / 3.));
    return curve;
}

std::vector<std::tuple<SurfacePoint, double>>
getCurveDistances(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
                  const std::vector<SurfacePoint>& curve) {
    std::vector<std::tuple<SurfacePoint, double>> distances{
        std::make_tuple(curve[0], 0)};

    double totalDist = 0;

    for (size_t iP = 0; iP + 1 < curve.size(); iP++) {
        SurfacePoint pPrev = curve[iP];
        SurfacePoint pCurr = curve[iP + 1];
        double dist =
            exactgeodesic::compute_surface_distance(geom, pPrev, pCurr);
        totalDist += dist;
        distances.push_back(std::make_tuple(pCurr, totalDist));
    }
    return distances;
}

// result is positive for vertices to the left of the curve and negative to the
// right
VertexData<double> diffuseLeft(ManifoldSurfaceMesh& mesh,
                               VertexPositionGeometry& geom,
                               VectorHeatMethodSolver& vhm,
                               const std::vector<SurfacePoint>& curve) {
    VertexData<double> source(mesh, 0);

    // TODO: handle curve start and end
    for (size_t iP = 1; iP + 2 < curve.size(); iP++) {
        SurfacePoint pPrev = curve[iP];
        SurfacePoint pNext = curve[iP + 1];
        Edge ePrev         = pPrev.edge;
        Edge eNext         = pNext.edge;

        Face f         = sharedFace(pPrev, pNext);
        Halfedge heOpp = f.halfedge();
        while (heOpp.edge() == ePrev || heOpp.edge() == eNext) {
            heOpp = heOpp.next();
        }
        // f's edges should now be ePrev, eNext, and heOpp.edge()

        double orientation = 1;
        if (heOpp.next().edge() == ePrev) {
            orientation = -1;
        }
        // orientation = 1 if heOpp is parallel to the curve, and -1 if it is
        // antiparallel


        double dist =
            exactgeodesic::compute_surface_distance(geom, pPrev, pNext);

        source[heOpp.next().tipVertex()] = orientation * dist;
        source[heOpp.tailVertex()]       = -orientation * dist;
        source[heOpp.tipVertex()]        = -orientation * dist;
    }

    return vhm.scalarDiffuse(source);
}

VertexData<std::pair<double, Vector2>>
curveLogMap(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
            const std::vector<SurfacePoint>& curve) {

    VectorHeatMethodSolver vhm(geom);
    std::vector<std::tuple<SurfacePoint, double>> zippedDistances =
        getCurveDistances(mesh, geom, curve);
    double curveLen                 = std::get<1>(zippedDistances.back());
    VertexData<double> closestPoint = vhm.extendScalar(zippedDistances);
    VertexData<double> leftCoeff    = diffuseLeft(mesh, geom, vhm, curve);

    // TODO: this factors matrices again, which is unnecessary. Use vhm instead?
    HeatMethodDistanceSolver hm(geom);
    VertexData<double> curveDist = hm.computeDistanceCurve(curve);
    // VertexData<double> oldDist   = hm.computeDistance(curve);
    // psMesh->addVertexDistanceQuantity("oldDist", oldDist);

    VertexData<Vector2> logStart = vhm.computeLogMap(curve.front());
    VertexData<Vector2> logEnd   = vhm.computeLogMap(curve.back());

    // Rotate vectors in log map vectors to be in basis of curve
    // (the basis where the tangent vector to the curve points in the x
    // direction)

    // TODO: currently assumes that the curve starts inside a face. Make more
    // general
    geom.requireHalfedgeVectorsInFace();
    auto startCoords = vertexCoordinatesInTriangle(geom, curve[0].face);
    Vector2 start    = baryCoordsToFaceCoords(startCoords, curve[0].faceCoords);
    Vector2 next     = baryCoordsToFaceCoords(
        startCoords, curve[1].inFace(curve[0].face).faceCoords);
    Vector2 startToCurve = (next - start).normalize().inv();
    int N                = curve.size() - 1;
    auto endCoords       = vertexCoordinatesInTriangle(geom, curve[N].face);
    Vector2 end  = baryCoordsToFaceCoords(endCoords, curve[N].faceCoords);
    Vector2 prev = baryCoordsToFaceCoords(
        endCoords, curve[N - 1].inFace(curve[N].face).faceCoords);
    Vector2 endToCurve = (end - prev).normalize().inv();
    geom.unrequireHalfedgeVectorsInFace();
    for (Vertex v : mesh.vertices()) {
        logStart[v] *= startToCurve;
        logEnd[v] *= endToCurve;
    }

    VertexData<std::pair<double, Vector2>> curveLog(mesh);
    VertexData<Vector2> justLog(mesh);
    VertexData<double> distToStart(mesh);
    VertexData<double> hmDistToStart = hm.computeDistance(curve.front());
    VertexData<double> hmDistToEnd   = hm.computeDistance(curve.back());
    VertexData<double> closestPointType(mesh);
    for (Vertex v : mesh.vertices()) {
        double dStart = hmDistToStart[v];
        double dEnd   = hmDistToEnd[v];
        double dMid   = curveDist[v];

        logStart[v] = logStart[v].normalize() * dStart;
        logEnd[v]   = logEnd[v].normalize() * dEnd;

        distToStart[v] = dStart;

        double tol = 1;
        if (dStart < dEnd && dStart <= dMid * tol) {
            curveLog[v]         = std::make_pair(0, logStart[v]);
            closestPointType[v] = -1;
        } else if (dEnd < dStart && dEnd <= dMid * tol) {
            curveLog[v]         = std::make_pair(1, logEnd[v]);
            closestPointType[v] = 1;
        } else {
            closestPointType[v] = 0;
            curveLog[v]         = std::make_pair(
                closestPoint[v] / curveLen,
                Vector2{0, copysign(curveDist[v], leftCoeff[v])});
        }
        justLog[v] = curveLog[v].second;
    }

    psMesh->addVertexScalarQuantity("closest point", closestPoint)
        ->setEnabled(true);
    // psMesh->addVertexScalarQuantity("onLeft", leftCoeff)->setEnabled(true);
    // psMesh->addVertexScalarQuantity("closest point type", closestPointType);
    psMesh->addVertexDistanceQuantity("dist to curve", curveDist)
        ->setEnabled(true);
    psMesh->addVertexDistanceQuantity("dist to start", distToStart);
    psMesh->addVertexDistanceQuantity("hm dist to start", hmDistToStart);
    // psMesh->addVertexParameterizationQuantity("log start", logStart)
    //     ->setStyle(polyscope::ParamVizStyle::LOCAL_RAD);
    geom.requireEdgeLengths();
    psMesh->addEdgeScalarQuantity("len", geom.edgeLengths);
    geom.unrequireEdgeLengths();
    psMesh->addVertexParameterizationQuantity("logs", justLog)
        ->setStyle(polyscope::ParamVizStyle::LOCAL_RAD)
        ->setEnabled(true);

    std::vector<double> plainDistances;
    for (const auto& zd : zippedDistances)
        plainDistances.push_back(std::get<1>(zd));
    psCurve->addNodeScalarQuantity("dist", plainDistances)->setEnabled(true);

    return curveLog;
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("Geometry program");
    args::Positional<std::string> inputFilename(parser, "mesh",
                                                "Mesh to be processed.");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string filename = "../../meshes/bunny_small.obj";
    // Make sure a mesh name was given
    if (inputFilename) {
        filename = args::get(inputFilename);
    }

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = myCallback;

    // Load mesh
    std::tie(mesh, geom) = readManifoldSurfaceMesh(filename);

    psMesh = polyscope::registerSurfaceMesh("mesh", geom->vertexPositions,
                                            mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));

    curve = generateCurve(*mesh, *geom);

    psCurve = drawCurve(*mesh, *geom, curve);

    curveLogMap(*mesh, *geom, curve);


    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
