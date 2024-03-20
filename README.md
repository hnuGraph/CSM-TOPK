# CSMTopK
## Introduction
 We creatively propose a new search process, called as NewSP, where items in the model are operations instead of extensions. We incorporate innovative operation-level expansion postponement in NewSP, which could avoid premature expansions without hindering the original pruning power of the selected matching order. With expansion postponement, there could be multiple consecutive expansions in NewSP model, and we take this unique opportu- nity to design a multi-expansion strategy for further optimization. Also, in our new model, candidate set may not be expanded imme- diately once determined, which provides another unique chance for us to design cache strategies for candidate set reuse. We additionally propose adaptive index filtering over our multi-expansion and candidate reuse for performance enhancement.

## Compile

Our framework requires c++17 and GCC 9.x (or later). Under the root directory of the project, execute the following commands to compile the source code.

```shell
mkdir build
cd build
cmake ..
make
```

## Execute

After a successful compilation, the binary file is created under the `build/` directory. One can execute CSMTOPK using the following command.

```shell
build/csm -q <query-graph-path> -d <data-graph-path> -u <update-stream-path> --ksize <the value of k>
```


### Commandline Parameters

Other commandline parameters supported by the framework are listed in the following table.

| Command Line Parameters | Description                                                 | Valid Value  | Default Value |
| ----------------------- | ----------------------------------------------------------- | ------------ | ------------- |
| --time-limit            | Time limit for the incremental matching phase (in seconds). | 0-4294967295 | 3600          |
| --report-initial        | Perform initial matching or not.                            | on/off       | on            |
| --initial-time-limit    | Time limit for the initial matching phase (in seconds).     | 0-4294967295 | 4294967295    |
| --print-prep            | Print preprocessing results or not.                         | on/off       | on  
--print-result|	print update top k  results or not|on/off|off|
| --print-enum            | print update top k  results or not                          | on/off       | off           |
| --ul                    | the size of deletion edges                                  | 0-10000      | 5000          |
| --qInfo                 | the path of Initial result of top k set                     |              | ""            |

For example, if one requires the framework (1)to print top k dense subgraphs results if $A^k_t$ changes.   (2) to spend at most 1 hour (3600 seconds) on the update stream, (3) $k=300$,then the command should be

```shell
build/csm -q <query-graph-path> -d <data-graph-path> -u <update-stream-path> --time-limit 3600 --print-result 1 --ksize 300
```

## Input File Format
Both the input query graph and data graph are vertex- and edge-labeled. Each edge on the data graph has a weight.  Each vertex is represented by a distinct unsigned integer (from 0 to 4294967295). There is at most one edge between two arbitrary vertices. 

### Query Graph

Each line in the query graph file represent a vertex or an edge.

1. A vertex is represented by `v <vertex-id> <vertex-label>`;
2. An edge is represented by `e <vertex-id-1> <vertex-id-2> <edge-label>`.

The two endpoints of an edge must appear before the edge. For example, 

```
v 0 0
v 1 0
v 2 1
e 0 1 0
e 0 2 1
e 2 1 2
```

### Initial Data Graph

Each line in the query graph file represent a vertex or an edge.

1. A vertex is represented by `v <vertex-id> <vertex-label>`;

2. An edge is represented by `e <vertex-id-1> <vertex-id-2> <edge-label> <edege-weight> `
```
v 0 0
v 1 0
v 2 1
e 0 1 0 1000
e 0 2 1 2000
e 2 1 2 3000
```

### Graph Update Stream

Graph update stream is a collection of insertions and deletions of a vertex or an edge.

1. A vertex insertion is represented by `v <vertex-id> <vertex-label>`;
2. A vertex deletion is represented by `-v <vertex-id> <vertex-label>`;
3. An edge insertion is represented by `e <vertex-id-1> <vertex-id-2> <edge-label> <edege-weight>` ;
4. An edge deletion is represented by `-e <vertex-id-1> <vertex-id-2> <edge-label> <edege-weight>`;

The vertex or edge to be deleted must exist in the graph, and the label must be the same as that in the graph. If an edge is inserted to the data graph, both its endpoints must exist. For example,

```
v 0 1
e 0 3 2 1000
-v 1 2
-e 0 1 0 1000
```

##  Datasets

We provide 4 datasets in our experiment 

1. Amazon dataset        [[link](https://snap.stanford.edu/data/com-Amazon.html)]
2. Livejournal dataset   [[link](https://snap.stanford.edu/data/soc-LiveJournal1.html)]
3. Human dataset       [[link](http://hprd.org/index_html)]
4. Youtube dataset       [[link](https://snap.stanford.edu/data/com-Youtube.html)]



**Summary of Datasets**

| **Datasets** |          **Type**           | **Vertexes** | **Edges**  | **L(V)** | **Average Degree** |
| :----------: | :-------------------------: | :----------: |:----------:|:--------:|:------------------:|
|    Amazon    |       Product network       |   403,394    | 1,015,000  |  6					  |   5.03       		    |
| Livejournal  |      Community network      |  4,847,571   | 30,005,000 |  30				  |       12.38        |
|    Human     | Protein interaction network |    4,674     |   81,282   |  44				  |       34.78        |
|   Youtube    |       Social network        |  1,134,890   | 2,015,000  |  25				  |        3.55        |