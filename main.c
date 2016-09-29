#include "header.h"

int read(Pixel*** v, FILE* in, FILE* out){
    /* Citeste matricea de pixeli din poza
     * Scrie headerul in output */
    char header[55];
    int size;
    int i, j;
    fread(header, 1, 54, in);
    fwrite(header, 1, 54, out);
    fseek(in, 18, SEEK_SET);
    fread(&size, 4, 1, in);
    fseek(in, 54, SEEK_SET);
    /* Alocare memorie pentru matrice de pixeli */
    *v = (Pixel **) malloc(size * sizeof(Pixel *));
    for (i = 0; i < size; i++)
        (*v)[i] = (Pixel *) malloc(size * sizeof(Pixel));
    for (i = size-1; i >= 0; i--)
        for (j = 0; j < size; j++)
            fread((*v)[i][j].colors, 1, 4, in);
    return size;
}

int readCompressed(QuadtreeNode**vector, Pixel ***v, FILE *in, FILE *out){
    /* Citirea datelor comprimate din fisier
     * Se scrie headerul in output
     * Se aloca memorie pentru vector si matrice
     * Se citeste vectorul*/

    char header[55];
    int size;
    int i, j;
    fread(header, 1, 54, in);
    fwrite(header, 1, 54, out);
    fseek(in, 18, SEEK_SET);
    fread(&size, 4, 1, in);
    fseek(in, 54, SEEK_SET);
    int Node, Leaf;
    fread(&Leaf, 4, 1, in);
    fread(&Node, 4, 1, in);
    *vector = (QuadtreeNode*) malloc(Node*sizeof(QuadtreeNode));
    fread(*vector, sizeof(QuadtreeNode), Node, in);
    *v = (Pixel **) malloc(size * sizeof(Pixel *));
    for (i = 0; i < size; i++)
        (*v)[i] = (Pixel *) malloc(size * sizeof(Pixel));
    return size;
}

int equals ( Pixel A, Pixel B){
    /* Verifica daca doua structuri de tip Pixel contin aceleasi valori */
    int i;
    for (i=0;i<=3;i++)
        if ( A.colors[i] != B.colors[i] )
            return 0;
    return 1;
}

int homogenous (Pixel** v,  int left,  int top,  int right,  int bottom) {
    /* Verifica daca patratul din regiunea left-right, top-bottom
     * contine aceeasi culoare */
    Pixel A = v[top][left];
    int i, j;
    for (i=top; i<=bottom; i++)
        for (j=left; j<=right; j++)
            if ( equals(A, v[i][j]) != 1 )
                return 0;
    return 1;
}

void buildQuadtree (Pixel ** v, Tree* root,  int left,  int top,  int right,  int bot){
    /* Construire arbore cuaternar */
    if ( homogenous(v, left, top, right, bot) == 1 ){
        /* Creare frunza in arbore */
        root->top_left = NULL;
        root->top_right = NULL;
        root->bottom_left = NULL;
        root->bottom_right = NULL;
        root->blue = v[top][left].colors[0];
        root->green = v[top][left].colors[1];
        root->red = v[top][left].colors[2];
        root->reserved = v[top][left].colors[3];
        root->area = (uint32_t) (bot-top+1)*(right-left+1);

    }
    else {
        /* Nodul nu este frunza. Se creeaza copiii si se porneste recursiv
         * in fiecare dintre ei */
        int Xmid =  (left+right)/2;
        int Ymid = (top+bot)/2;
        root->area = (uint32_t) (bot - top + 1) * (right - left + 1);
        root->top_left = (Tree *) calloc(1, sizeof(Tree));
        root->top_right = (Tree *) calloc(1, sizeof(Tree));
        root->bottom_right = (Tree *) calloc(1, sizeof(Tree));
        root->bottom_left = (Tree *) calloc(1, sizeof(Tree));
        buildQuadtree(v, root->top_right, left, top, Xmid, Ymid);
        buildQuadtree(v, root->bottom_right, Xmid+1, top, right, Ymid);
        buildQuadtree(v, root->bottom_left, Xmid+1, Ymid+1, right, bot);
        buildQuadtree(v, root->top_left, left, Ymid + 1, Xmid, bot);
    }
}

void assignTreeToVector ( Tree root, QuadtreeNode* element){
    /* Copiere element din arbore in vector */
    element->area = root.area;
    element->green = root.green;
    element->blue = root.blue;
    element->red = root.red;
    element->reserved = root.reserved;
    /* Daca nodul este frunza, se scriu indicii -1 in dreptul copiilor */
    if (root.top_left == NULL){
        element->top_left = -1;
        element->top_right = -1;
        element->bottom_right = -1;
        element->bottom_left = -1;
    }
}

void assignVectorToTree(Tree *element, QuadtreeNode *root) {
    /* Copiere element din vector in arbore */
    element->area = root->area;
    element->green = root->green;
    element->blue = root->blue;
    element->red = root->red;
    element->reserved = root->reserved;
    if (root->top_left == -1) {
        element->top_left = NULL;
        element->top_right = NULL;
        element->bottom_right = NULL;
        element->bottom_left = NULL;
    }
}

int TreeToVector (Tree root, QuadtreeNode* vector, int pos){
    /* Transformarea arborelui in vector.
     * Se adauga recursiv mai intai patratul din stanga sus, apoi cele din
     * dreapta sus, dreapta jos si stanga jos. */
    int k = pos;
    assignTreeToVector(root, vector+pos);
    if (root.top_left!=NULL && root.top_right!=NULL && root.bottom_right !=NULL && root.bottom_left!=NULL){
        k=k+1;
        vector[pos].top_left =  k;
        k = TreeToVector(*root.top_left, vector, k);
        k = k + 1;
        vector[pos].top_right = k;
        k = TreeToVector(*root.top_right, vector, k);
        k = k + 1;
        vector[pos].bottom_right = k;
        k = TreeToVector(*root.bottom_right, vector, k);
        k = k + 1;
        vector[pos].bottom_left = k;
        k = TreeToVector(*root.bottom_left, vector, k);
    }
    return k;
}

void copyPixel (Pixel **v, int left, int top, int right, int bot, Pixel color){
    /* Copiaza pixelul transmis ca parametru in matrice intre pozitiile
     * left - right, top - bottom */
    int i, j;
    for (i=left; i<=right;i++)
        for (j=top;j<=bot;j++) {
            v[i][j].colors[0]=color.colors[0];
            v[i][j].colors[1] = color.colors[1];
            v[i][j].colors[2] = color.colors[2];
            v[i][j].colors[3] = color.colors[3];
         }
}

void buildMatrix (Pixel **v, Tree *root, int left, int top, int right, int bot){
    /* Construire matrice din arbore cuaternar */
    Pixel aux;
    if (root->top_left==NULL){
        /* Daca root este frunza, se scrie culoarea in matrice */
        aux.colors[0] = root->blue;
        aux.colors[1] = root->green;
        aux.colors[2] = root->red;
        aux.colors[3] = root->reserved;
        copyPixel(v,left,top,right,bot,aux);
    }
    else {
        /* Altfel, se pleaca recursiv in cei 4 fii pana la gasirea frunzelor */
        int Xmid = (left + right) / 2;
        int Ymid = (top + bot) / 2;
        buildMatrix(v, root->top_left, left, top, Xmid, Ymid);
        buildMatrix(v, root->top_right, Xmid + 1, top, right, Ymid);
        buildMatrix(v, root->bottom_right, Xmid + 1, Ymid + 1, right, bot);
        buildMatrix(v, root->bottom_left, left, Ymid + 1, Xmid, bot);
    }

}
void VectorToTree ( Tree ** root, QuadtreeNode* vector, int pos){
    /* Construieste arborele cuaternar plecand de la vector */
    *root = (Tree*) malloc (sizeof(Tree));
    /* Se copiaza valorile din vector in root */
    assignVectorToTree(*root,vector+pos);
    if (vector[pos].top_left != -1) {
        /* Se construiesc recursiv nodurile corespunzatoare celor 4 fii */
       VectorToTree(&((*root)->top_left), vector, vector[pos].top_left);
       VectorToTree(&((*root)->top_right), vector, vector[pos].top_right);
       VectorToTree(&((*root)->bottom_right), vector, vector[pos].bottom_right);
       VectorToTree(&((*root)->bottom_left), vector, vector[pos].bottom_left);
    }
}
void rotire (QuadtreeNode * element, int k){
    /* Rotirea pozei. Se permuta indicii fiilor pentru fiecare element din vector*/
    int aux;
    if (element[k].top_left != -1 ){
        aux = element[k].top_left;
        element[k].top_left= element[k].top_right;
        element[k].top_right= element[k].bottom_right;
        element[k].bottom_right= element[k].bottom_left;
        element[k].bottom_left=aux;
        rotire(element, element[k].top_left);
        rotire(element, element[k].top_right);
        rotire(element, element[k].bottom_right);
        rotire(element, element[k].bottom_left);
    }
}
void eliberare(Tree * root){
    /* Eliberare memorie arbore */
    if (root == NULL)
        return;
    else {
        eliberare(root->top_left);
        eliberare(root->top_right);
        eliberare(root->bottom_left);
        eliberare(root->bottom_right);
        free(root);
    }
}
int isAncestor ( Tree * Node, Pixel A){
    /* Verifica daca nodul este stramos al unui pixel A transmis ca parametru */
    if (Node->top_left == NULL){
        if (Node->blue == A.colors[0] && Node->green == A.colors[1] && Node->red == A.colors[2] )
            return 1;
        else return 0;
    }
    /* Daca nodul nu este frunza de culoarea pixelului A, se cauta recursiv
     * in cei 4 fii */
    return isAncestor(Node->top_left,A) + isAncestor(Node->top_right,A) + isAncestor(Node->bottom_left, A) + isAncestor(Node->bottom_right, A);
}

void LCA ( Tree * Node , Pixel A, Pixel B, Tree** Anc){
    /* Se cauta cel mai mic stramos al pixelilor A si B
     * Daca unul dintre fii este stramos al celor doua culori, se pleaca recursiv
     * pe acea directie in cautarea celui mai mic stramos comun.
     * De fiecare data cand se gaseste un stramos cu o suprafata mai mica
     * decat a celui curent, se scrie in Anc */

    if (isAncestor(Node->top_left, A) != 0 && isAncestor(Node->top_left, B) != 0) {
               if (Node->top_left->area < (*Anc)->area)
                   *Anc = Node->top_left;
               LCA (Node->top_left, A,B, Anc);
           }
    if (isAncestor(Node->top_right, A) != 0 && isAncestor(Node->top_right, B) != 0) {
        if (Node->top_right->area < (*Anc)->area)
            *Anc = Node->top_right;
        LCA(Node->top_right, A, B, Anc);
    }
    if (isAncestor(Node->bottom_left, A) != 0 && isAncestor(Node->bottom_left, B) != 0) {
        if (Node->bottom_left->area < (*Anc)->area)
            *Anc = Node->bottom_left;
        LCA(Node->bottom_left, A, B, Anc);
    }
    if (isAncestor(Node->bottom_right, A) != 0 && isAncestor(Node->bottom_right, B) != 0) {
        if (Node->bottom_right->area < (*Anc)->area)
            *Anc = Node->bottom_right;
        LCA(Node->bottom_right, A, B, Anc);

    }
}

uint32_t count_nodes (Tree *root){
    if (root->top_left == NULL)
        return 1;
    else return 1 + count_nodes(root->top_left) + count_nodes(root->top_right) + count_nodes(root->bottom_left) + count_nodes(root->bottom_right);
}

uint32_t count_leaves(Tree * root){
    if (root->top_left == NULL)
        return 1;
    else
        return count_leaves(root->top_left) + count_leaves(root->top_right) + count_leaves(root->bottom_left) + count_leaves(root->bottom_right);
}
int main (int argc,  char* argv[]){
    FILE * in, *out;
    Pixel** Matrix;
    QuadtreeNode *vector;
    int size;
    int old_size;
    int i, j;
    uint32_t Node = 0, Leaf = 0;
    /* Compresie imagine */
    if (strcmp(argv[1], "-c") == 0) {
        in = fopen(argv[2], "rb");
        out = fopen(argv[3], "wb");
        size = read(&Matrix, in, out);
        old_size = size;
        Tree *root = (Tree *) calloc(1,sizeof(Tree));
        buildQuadtree(Matrix, root, 0, 0, size - 1, size - 1);
        Node = count_nodes(root);
        Leaf = count_leaves(root);
        vector = (QuadtreeNode *) calloc(Node, sizeof(QuadtreeNode));
        TreeToVector(*root, vector, 0);
        fwrite(&Leaf, 4, 1, out);
        fwrite(&Node, 4, 1, out);
        fwrite(vector, sizeof(QuadtreeNode), Node, out);
        eliberare(root);
    }
    /* Decompresie imagine */

    else if (strcmp(argv[1], "-d") == 0) {
        in = fopen(argv[2], "rb");
        out = fopen(argv[3], "wb");
        size = readCompressed(&vector, &Matrix, in, out);
        old_size=size;
        Tree *root;
        VectorToTree(&root, vector, 0);
        buildMatrix(Matrix, root, 0, 0, size - 1, size - 1);
        for (i = 0; i < size; i++)
            for (j = 0; j < size; j++)
                fwrite(Matrix[i][j].colors, 1, 4, out);
        eliberare(root);
    }
    /* Rotire imagine */
    else if (strcmp(argv[1], "-r") == 0) {
        in = fopen(argv[3], "rb");
        out = fopen(argv[4], "wb");
        int nr_rot;
        sscanf(argv[2],"%d", &nr_rot);
        size = read(&Matrix, in, out);
        old_size = size;
        Tree *root = (Tree *) malloc(sizeof(Tree));
        buildQuadtree(Matrix, root, 0, 0, size - 1, size - 1);
        Node = count_nodes(root);
        Leaf = count_leaves(root);
        vector = (QuadtreeNode *) malloc(Node * sizeof(QuadtreeNode));
        TreeToVector(*root, vector, 0);
        //eliberare(root);
        for (i=1;i<=nr_rot;i++){
            rotire(vector, 0);
        }
        VectorToTree(&root, vector, 0);
        buildMatrix(Matrix, root, 0, 0, size - 1, size - 1);
        for (i = 0; i < size; i++)
            for (j = 0; j < size; j++)
                fwrite(Matrix[i][j].colors, 1, 4, out);
        //eliberare(root);
    }
    /* Bonus */
    else if (strcmp(argv[1], "-b") == 0){
        int r1,b1,g1,r2,b2,g2;
        sscanf(argv[2], "%d", &r1);
        sscanf(argv[3], "%d", &g1);
        sscanf(argv[4], "%d", &b1);
        sscanf(argv[5], "%d", &r2);
        sscanf(argv[6], "%d", &g2);
        sscanf(argv[7], "%d", &b2);
        in = fopen(argv[8], "rb");
        out = fopen(argv[9], "wb");
        size = read(&Matrix, in, out);
        old_size =  size;
        Tree *root = (Tree *) malloc(sizeof(Tree));
        buildQuadtree(Matrix, root, 0, 0, size - 1, size - 1);
        Node = count_nodes(root);
        Leaf = count_leaves(root);
        Pixel A, B;
        A.colors[0] = (unsigned char) b1;
        A.colors[1] = (unsigned char) g1;
        A.colors[2] = (unsigned char) r1;
        B.colors[0] = (unsigned char) b2;
        B.colors[1] = (unsigned char) g2;
        B.colors[2] = (unsigned char) r2;
        Tree * Stramos = root;
        LCA ( root, A, B, &Stramos);
        vector = (QuadtreeNode *) malloc(Node * sizeof(QuadtreeNode));
        TreeToVector(*Stramos, vector, 0);
        eliberare(root);
        VectorToTree(&Stramos, vector, 0);
        //eliberare(Stramos);
        size = 1;
        while(size*size != vector[0].area){
            size=size*2;
        }
        fseek(out, 18, SEEK_SET);
        int new_size = (int) size;
        fwrite(&new_size, 4, 1, out);
        fwrite(&new_size, 4, 1, out);
        fseek(out, 54, SEEK_SET);
        buildMatrix(Matrix, Stramos, 0, 0, size - 1, size - 1);
        for (i = 0; i < size; i++)
            for (j = 0; j < size; j++)
                fwrite(Matrix[i][j].colors, 1, 4, out);

    }
    for (i = 0; i < old_size; i++)
        free(Matrix[i]);
        free(Matrix);
        free(vector);
        fclose(in);
        fclose(out);

}

