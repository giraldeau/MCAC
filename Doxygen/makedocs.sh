#!/bin/sh

FLOWGEN_DIR="../../Flowgen/"
SRCS="../Generateur_DLCA/ ../Analyse_DLCA ../Analyse_Projections2D"


function flowgen
{
    SRC_DIR="$1"
    SRC_NAME=$(basename $SRC_DIR)


    #List sources files
    source_files=$(ls ${SRC_DIR}/*.cpp)

    #Make flowgen documentation
    for file in ${source_files}; do
        python3 ${FLOWGEN_DIR}/build_db.py ${file};
    done
    for file in ${source_files}; do
        python3 ${FLOWGEN_DIR}/makeflows.py ${file};
    done
    java -jar ${FLOWGEN_DIR}/plantuml.jar flowdoc/aux_files/*.txt
    #cat <<EOF > flowdoc/aux_files/runphase
    cp -r ${FLOWGEN_DIR}/htmlCSSandJS/ flowdoc/htmlCSSandJS/
    for file in ${source_files}; do
        python3 ${FLOWGEN_DIR}/makehtml.py ${file};
    done
    if [ -d $SRC_NAME"_flowdoc" ] ; then rm -r $SRC_NAME"_flowdoc";fi
    mv flowdoc $SRC_NAME"_flowdoc"
}


for SRC_DIR in $SRCS;do
    flowgen $SRC_DIR
done
#Make Doxygen documentation
doxygen Doxyfile
