package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2015 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */


import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;


/**
 *
 * @author sc
 */
public class RnaseqExpressionConverter extends BioFileConverter
{
    //
    private static final String DATASET_TITLE = "RNAseq expression";
    private static final String DATA_SOURCE_NAME = "whoknows";
    private static final Logger LOG = Logger.getLogger(RnaseqExpressionConverter.class);
    private static final int STUDIES_NR = 4;
    private Item org;
    // TODO get taxid from property file
    private static final String TAX_ID = "3702";

    private static Map<String, String> studies = null;
    private Map<String, String> geneItems = new HashMap<String, String>();
    private Map<String, String> transcriptItems = new HashMap<String, String>();

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public RnaseqExpressionConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        try {
            createOrganismItem();
        } catch (ObjectStoreException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    /**
     *
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
        // There are two types of files:
        // The following code works out which file we are reading
        File currentFile = getCurrentFile();

        if ("4samples.gene".equals(currentFile.getName())) {
            LOG.info("FILE 4 GENES");
            processFile(reader, "gene", org);
        } else if ("4samples.transcript".equals(currentFile.getName())) {
            LOG.info("FILE 4 TRANSCRIPTS");
            processFile(reader, "transcript", org);
        } else {
            throw new IllegalArgumentException("Unexpected file: "
                    + currentFile.getName());
        }
    }

    /**
     * Process all rows of the expression files
     *
     * @param reader
     *            a reader for the
     *            SAMPLE_NRsample.gene and SAMPLE_NRsample.transcript
     *            files
     * @throws IOException
     * @throws ObjectStoreException
     */
    private void processFile(Reader reader, String type, Item organism)
        throws IOException, ObjectStoreException {
        Iterator<?> tsvIter;
        try {
            tsvIter = FormattedTextParser.parseTabDelimitedReader(reader);
        } catch (Exception e) {
            throw new BuildException("cannot parse file: " + getCurrentFile(), e);
        }
        String [] headers = null;
        int lineNumber = 0;

        studies = new HashMap<String, String>();

        while (tsvIter.hasNext()) {
            String[] line = (String[]) tsvIter.next();
            LOG.info("SCOREg " + line[0]);

            if (lineNumber == 0) {
                // column headers - strip off any extra columns
                int end = 0;
                for (int i = 0; i < line.length; i++) {
                    if (StringUtils.isEmpty(line[i])) {
                        break;
                    }
                    end++;
                }
                headers = new String[end];
                System.arraycopy(line, 0, headers, 0, end);
                LOG.info("HHH " + headers.length + " end=" + end);
            } else {
                String primaryId = line[0]; //Gene id
                // if empty lines at the end of the file
                if (StringUtils.isEmpty(primaryId)) {
                    break;
                }

                if (type.equalsIgnoreCase("gene")) {
                    createBioEntity(primaryId, "Gene");
                }
                if (type.equalsIgnoreCase("transcript")) {
                    createBioEntity(primaryId, "Transcript");
                }

                // scores start from column 2 and end at STUDIES_NR +1 which is headers[1,SN]
                for (int i = 1; i <= STUDIES_NR; i++) {
                    String col = headers[i].replace("_TMP", "");
                    //                    col = correctOfficialName(col);
                    if (!studies.containsKey(col)) {
                        Item experiment = createExperiment(col);
                        studies.put(col, experiment.getIdentifier());
                    }
                    Item score = createRNASeqExpression(line[i], type);
                    if (type.equalsIgnoreCase("gene")) {
                        score.setReference("gene", geneItems.get(primaryId));
                    }
                    if (type.equalsIgnoreCase("transcript")) {
                        score.setReference("transcript", transcriptItems.get(primaryId));
                    }
                    score.setReference("study", studies.get(col));
                    score.setReference("organism", organism);
                    store(score);
                }
            }
            lineNumber++;
        }
    }

    /**
     * Create and store a RnaseqExpression item on the first time called.
     *
     * @param score the expression score
     * @param type gene or transcript
     * @return an Item representing the GeneExpressionScore
     */
    private Item createRNASeqExpression(String score, String type) throws ObjectStoreException {
        Item expression = createItem("RnaseqExpression");
        expression.setAttribute("TPM", score);
        expression.setAttribute("type", type);
        return expression;
    }

    /**
     * Create and store a BioEntity item on the first time called.
     *
     * @param primaryId the primaryIdentifier
     * @param type gene or transcript
     * @throws ObjectStoreException
     */
    private void createBioEntity(String primaryId, String type) throws ObjectStoreException {
        Item bioentity = null;
        LOG.info("BIO: " + type + " -- " + primaryId);
        if ("Gene".equals(type)) {
            if (!geneItems.containsKey(primaryId)) {
                bioentity = createItem("Gene");
                bioentity.setAttribute("primaryIdentifier", primaryId);
                store(bioentity);
                geneItems.put(primaryId, bioentity.getIdentifier());
            }
        } else if ("Transcript".equals(type)) {
            if (!transcriptItems.containsKey(primaryId)) {
                bioentity = createItem("Transcript");
                bioentity.setAttribute("primaryIdentifier", primaryId);
                store(bioentity);
                transcriptItems.put(primaryId, bioentity.getIdentifier());
            }
        }
    }

    /**
     * Create and store a organism item on the first time called.
     *
     * @throws ObjectStoreException os
     */
    protected void createOrganismItem() throws ObjectStoreException {
        org = createItem("Organism");
        org.setAttribute("taxonId", TAX_ID);
        store(org);
    }

    /**
     * Create and store an Experiment item on the first time called.
     *
     * @param name the cell line name
     * @return an Item representing the CellLine
     */
    private Item createExperiment(String name) throws ObjectStoreException {
        Item e = createItem("Experiment");
        e.setAttribute("title", name);
        store(e);
        return e;
    }

}

