#include "mainwindow.h"
#include "ui_mainwindow.h"


/* **** début de la partie à compléter **** */
float MainWindow::faceArea(MyMesh* _mesh, int faceID)
{
    FaceHandle fh = _mesh->face_handle ( faceID );

    std::vector<VertexHandle> vertexes;

    MyMesh::FaceVertexIter fh_v = _mesh->fv_iter(fh);
    for(; fh_v.is_valid(); ++fh_v)
        vertexes.push_back ( *fh_v );

    float valuesAB[3];
    valuesAB[0] = _mesh->point(vertexes[1])[0] - _mesh->point(vertexes[0])[0];
    valuesAB[1] = _mesh->point(vertexes[1])[1] - _mesh->point(vertexes[0])[1];
    valuesAB[2] = _mesh->point(vertexes[1])[2] - _mesh->point(vertexes[0])[2];

    float valuesAC[3];
    valuesAC[0] = _mesh->point(vertexes[2])[0] - _mesh->point(vertexes[0])[0];
    valuesAC[1] = _mesh->point(vertexes[2])[1] - _mesh->point(vertexes[0])[1];
    valuesAC[2] = _mesh->point(vertexes[2])[2] - _mesh->point(vertexes[0])[2];

    VectorT<float, 3> vectorAB(valuesAB);
    VectorT<float, 3> vectorAC(valuesAC);

    VectorT<float, 3> product = vectorAB * vectorAC;
    float norm = product.norm();

    return norm / 2.0f;
}

float MainWindow::baryArea(MyMesh* _mesh, int vertID){
    float baryArea = 0;

    VertexHandle vh = _mesh->vertex_handle ( vertID );
    MyMesh::VertexFaceIter vf = _mesh->vf_iter ( vh );
    for ( ; vf.is_valid ( ) ; ++vf ) {
        FaceHandle current = *vf;
        baryArea += faceArea ( _mesh , current.idx( ) );
    }
    return baryArea / 3.0f;
}

float MainWindow::angleFF(MyMesh* _mesh, int faceID0,  int faceID1, int vertID0, int vertID1)
{
    int sign = 0;

    VertexHandle vh0 = _mesh->vertex_handle ( vertID0 );
    FaceHandle fh0 = _mesh->face_handle ( faceID0 );

    MyMesh::FaceVertexCWIter fh_cwv = _mesh->fv_cwiter ( fh0 );
    while ( fh_cwv.is_valid ( ) && *fh_cwv != vh0 ) ++fh_cwv;

    VertexHandle next = *++fh_cwv;

    if ( next.idx ( ) == vertID1 ) sign = -1;
    else sign = 1;

    OpenMesh::Vec3f normal0 (_mesh->normal ( fh0 ) );
    OpenMesh::Vec3f normal1 (_mesh->normal ( _mesh->face_handle ( faceID1 ) ) );

    float scalar = normal0 | normal1;

    return sign * acos ( scalar );
}

float MainWindow::angleEE(MyMesh* _mesh, int vertexID,  int faceID)
{
    FaceHandle fh = _mesh->face_handle ( faceID );
    VertexHandle vh = _mesh->vertex_handle ( vertexID );
    std::vector<VertexHandle> vertexes;

    MyMesh::FaceVertexIter fh_v = _mesh->fv_iter(fh);
    for(; fh_v.is_valid(); ++fh_v) {
        VertexHandle current = *fh_v;
        if( current.idx() != vertexID )
            vertexes.push_back ( current );
    }

    float valuesAB[3];
    valuesAB[0] = _mesh->point(vertexes[0])[0] - _mesh->point(vh)[0];
    valuesAB[1] = _mesh->point(vertexes[0])[1] - _mesh->point(vh)[1];
    valuesAB[2] = _mesh->point(vertexes[0])[2] - _mesh->point(vh)[2];

    float valuesAC[3];
    valuesAC[0] = _mesh->point(vertexes[1])[0] - _mesh->point(vh)[0];
    valuesAC[1] = _mesh->point(vertexes[1])[1] - _mesh->point(vh)[1];
    valuesAC[2] = _mesh->point(vertexes[1])[2] - _mesh->point(vh)[2];

    VectorT<float, 3> normalizedAB(VectorT<float, 3>(valuesAB).normalize());
    VectorT<float, 3> normalizedAC(VectorT<float, 3>(valuesAC).normalize());

    return acos ( normalizedAB | normalizedAC );
}

void MainWindow::H_Curv(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++) {
        MyMesh::VertexHandle current = *curVert;
        float val = 0.0f;

        for (MyMesh::VertexEdgeIter currentEdge = _mesh->ve_iter ( current ); currentEdge.is_valid(); currentEdge++)
        {
            MyMesh::EdgeHandle eh = *currentEdge;
            MyMesh::HalfedgeHandle heh0 = _mesh->halfedge_handle(eh, 0);
            MyMesh::HalfedgeHandle heh1 = _mesh->halfedge_handle(eh, 1);

            FaceHandle fh0 = _mesh->face_handle(heh0);
            FaceHandle fh1 = _mesh->face_handle(heh1);

            // Si l'arête est en bordure, on ne traite qu'une face
            if ( fh1.idx ( ) > _mesh->n_faces ( ) )
                fh1 = fh0;

            // Détermine l'autre sommet
            int vertex2ID = _mesh->to_vertex_handle(heh1).idx();
            if (vertex2ID == current.idx ( ) )
                vertex2ID = _mesh->to_vertex_handle(heh0).idx();

            // ||e_ij||
            //VectorT<float,3> e = _mesh->point(_mesh->vertex_handle(vertex2ID)) - _mesh->point(vh);
            OpenMesh::Vec3f currentOppVector = _mesh->point ( _mesh->vertex_handle ( vertex2ID ) ) - _mesh->point ( current );

            OpenMesh::Vec3f normal0 ( _mesh->normal ( fh0 ) );
            OpenMesh::Vec3f normal1 ( _mesh->normal ( fh1 ) );

            if ( ( ( normal0 % normal1 ) | currentOppVector ) < 0 )
            {
                FaceHandle tempF = fh0;
                fh0 = fh1;
                fh1 = tempF;
            }

            val += currentOppVector.norm ( ) * angleFF ( _mesh , fh0.idx ( ) , fh1.idx ( ) , current.idx() , vertex2ID );
        }
        val /= ( 4 * baryArea ( _mesh , current.idx ( ) ) );
        _mesh->data ( current ).value = val;
    }
}

void MainWindow::K_Curv(MyMesh* _mesh)
{
    for ( MyMesh::VertexIter curVert = _mesh->vertices_begin() ; curVert!=_mesh->vertices_end() ; ++curVert ) {
        VertexHandle current = *curVert;

        float area = baryArea ( _mesh , current.idx() );
        float angleEESum = 0;

        MyMesh::VertexFaceIter vf = _mesh->vf_iter ( current );
        for ( ; vf.is_valid ( ) ; ++vf ) {
            FaceHandle currentFace = *vf;
            angleEESum += angleEE ( _mesh , current.idx ( ) , currentFace.idx ( ) );
        }

        _mesh->data ( current ).value = ( 1 / area ) * ( 2 * M_PI - angleEESum );
    }
}


/* **** début de la partie boutons et IHM **** */
void MainWindow::on_pushButton_H_clicked()
{
    H_Curv(&mesh);
    displayMesh(&mesh, true); // true permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}

void MainWindow::on_pushButton_K_clicked()
{
    K_Curv(&mesh);
    displayMesh(&mesh, true); // true permet de passer en mode "carte de temperatures", avec une gestion automatique de la couleur (voir exemple)
}
/*
    Cette fonction est à utiliser UNIQUEMENT avec le fichier testAngleArea.obj
    Elle est appelée par le bouton "Test angles/aires"

    Elle permet de vérifier les fonctions faceArea, angleFF et angleEE.
    Elle doit afficher :

    Aire de la face 0 : 2
    Aire de la face 1 : 2
    Angle entre les faces 0 et 1 : 1.5708
    Angle entre les faces 1 et 0 : -1.5708
    Angle au sommet 1 sur la face 0 : 0.785398
*/

void MainWindow::on_pushButton_angleArea_clicked()
{
    qDebug() << "Aire de la face 0 :" << faceArea(&mesh, 0);
    qDebug() << "Aire de la face 1 :" << faceArea(&mesh, 1);

    qDebug() << "Aire barycentrique sous le sommet 1: " << baryArea ( &mesh , 1 );

    qDebug() << "Angle entre les faces 0 et 1 :" << angleFF(&mesh, 0, 1, 1, 2);
    qDebug() << "Angle entre les faces 1 et 0 :" << angleFF(&mesh, 1, 0, 1, 2);

    qDebug() << "Angle au sommet 1 sur la face 0 :" << angleEE(&mesh, 1, 0);
    qDebug() << "Angle au sommet 3 sur la face 1 :" << angleEE(&mesh, 3, 1);
}

void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}
/* **** fin de la partie boutons et IHM **** */

/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    vertexSelection = -1;
    edgeSelection = -1;
    faceSelection = -1;

    modevoisinage = false;

    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}
