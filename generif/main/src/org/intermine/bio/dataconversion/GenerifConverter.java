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

import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

/**
 *
 * @author
 */
public class GenerifConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(GenerifConverter.class);

    private static final String DATASET_TITLE = "GeneRIF";
    private static final String DATA_SOURCE_NAME = "NCBI";

    //
    private Item org;

    private static final String ATH_TAXID = "3702";

    private Map<String, String> pubItems = new HashMap<String, String>();
    private Map<String, String> geneItems = new HashMap<String, String>();

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     */
    public GenerifConverter(ItemWriter writer, Model model)
            throws ObjectStoreException {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
        createOrganismItem();
    }

    /**
     *
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
        File currentFile = getCurrentFile();

        if ("generifs_basic".equals(currentFile.getName())) {
            processScoreFile(reader, org);
        } else {
            LOG.info("WWSS skipping file: " + currentFile.getName());
            //            throw new IllegalArgumentException("Unexpected file: "
            //          + currentFile.getName());
        }
    }

    /**
     * Process all rows of the
     * Drosophila_Cell_Lines_and_Developmental_Stages_Gene_Scores.txt file
     *
     * @param reader
     *            a reader for the
     *            Drosophila_Cell_Lines_and_Developmental_Stages_Gene_Scores.txt
     *            file
     *
     * @throws IOException
     * @throws ObjectStoreException
     */
    private void processScoreFile(Reader reader, Item organism)
            throws IOException, ObjectStoreException {
        Iterator<?> tsvIter;
        try {
            tsvIter = FormattedTextParser.parseTabDelimitedReader(reader);
        } catch (Exception e) {
            throw new BuildException("cannot parse file: " + getCurrentFile(), e);
        }
        IdResolver athResolver = IdResolverService.getIdResolverByOrganism("3702");
        String pid = null;


        String [] headers = null;
        int lineNumber = 0;
        LOG.info("WW  BEGIN -----------------------");

        while (tsvIter.hasNext()) {
            String[] line = (String[]) tsvIter.next();

            //            if (lineNumber == 0) {
            //LOG.info("SCOREg " + line[0]);
            // column headers - strip off any extra columns - FlyAtlas
            // not necessary for expressionScore, but OK to keep the code
            //                int end = 0;
            //                for (int i = 0; i < line.length; i++) {
            //                    if (StringUtils.isEmpty(line[i])) {
            //                        break;
            //                    }
            //                    end++;
            //                }
            //                headers = new String[end];
            //                System.arraycopy(line, 0, headers, 0, end);
            //                LOG.info("WW header lenght " + headers.length);
            //                lineNumber++;
            //                continue;
            //            }

            String taxid = line[0];
            String geneId = line [1];
            String pubMedId = line [2];
            String timeStamp = line[3];
            String annotation = line [4];

            if (!taxid.equalsIgnoreCase(ATH_TAXID)) {
                continue;
            }

            int resCount = athResolver.countResolutions(taxid, geneId);
            if (resCount != 1) {
                LOG.info("RESOLVER: failed to resolve gene to one identifier, ignoring gene: "
                        + geneId + " count: " + resCount);
                continue;
            }

            // NOT WORKING!
            //          if (pid == null) {
            //              LOG.info("MISSING ID: " + geneId);
            //              continue;
            //          }


            pid = athResolver.resolveId(taxid, geneId).iterator().next();
            LOG.info("READING " + taxid + ": " + pid + "<->" + geneId + "|" + pubMedId + "|"
                    + timeStamp + "--" + annotation);

            Item ann = createGeneRIF(annotation, timeStamp);
            createBioEntity(pid, "Gene");
            createPublication(pubMedId);
            ann.setReference("gene", geneItems.get(pid));
            ann.setReference("organism", organism);
            ann.setReference("publication", pubItems.get(pubMedId));
            store(ann);

            lineNumber++;
        }
    }


    /**
     * Unify variations on similar official names.
     *
     * TODO, data from marc
     *
     * @param name the original 'official name' value
     * @param type cell line or developmental stage
     * @return a unified official name
     */
    private String correctOfficialName(String name, String type) {
        if (name == null) {
            return null;
        }

        if (type.equals(ATH_TAXID)) {
            name = name.replace("_", " ");

            if (name.matches("^emb.*\\d-\\dh$")) {
                name = name.replaceFirst("emb", "Embryo");
                name = name.replaceFirst("h", " h");
            }
            // Assume string like "L3_larvae_dark_blue" has the official name
            // "L3 stage larvae dark blue"
            if (name.matches("^L\\d.*larvae.*$")) {
                name = name.replace("larvae", "stage larvae");
            }
            // TODO "WPP_2days" is not in the database
            if (name.matches("^WPP.*$")) {
                if (name.endsWith("hr")) {
                    String[] strs = name.split(" ");
                    StringBuffer sb = new StringBuffer();
                    sb.append(strs[0]).append(" + ").append(strs[1]);
                    name = name.replaceFirst("hr", " h");
                } else if (name.endsWith("days")) {

                }
                name = name.replaceFirst("WPP", "White prepupae (WPP)");
            }
        }

        return name;
    }



    /**
     * Create and store a GeneRIF item on the first time called.
     *
     * @param annotation the RIF note
     * @param timeStamp
     * @return an Item representing the geneRIF
     */
    private Item createGeneRIF(String annotation, String timeStamp) throws ObjectStoreException {
        Item generif = createItem("Generif");
        generif.setAttribute("annotation", annotation);
        generif.setAttribute("timeStamp", timeStamp);

        return generif;
    }


    /**
     * Create and store a BioEntity item on the first time called.
     *
     * @param primaryId the primaryIdentifier
     * @param type gene or exon
     * @throws ObjectStoreException
     */
    private void createBioEntity(String primaryId, String type) throws ObjectStoreException {
        Item bioentity = null;

        if ("Gene".equals(type)) {
            if (!geneItems.containsKey(primaryId)) {
                bioentity = createItem("Gene");
                bioentity.setAttribute("symbol", primaryId);
                store(bioentity);
                geneItems.put(primaryId, bioentity.getIdentifier());
            }
        }
    }
    /**
     * Create and store a BioEntity item on the first time called.
     *
     * @param primaryId the primaryIdentifier
     * @param type gene or exon
     * @throws ObjectStoreException
     */
    private void createPublication(String primaryId) throws ObjectStoreException {
        Item pub = null;
        if (!pubItems.containsKey(primaryId)) {
            pub = createItem("Publication");
            pub.setAttribute("pubMedId", primaryId);
            store(pub);
            pubItems.put(primaryId, pub.getIdentifier());
        }
    }

//    private String getPub(String pubMedId) throws ObjectStoreException {
//        String itemId = pubs.get(pubMedId);
//        if (itemId == null) {
//            Item pub = createItem("Publication");
//            pub.setAttribute("pubMedId", pubMedId);
//            itemId = pub.getIdentifier();
//            pubs.put(pubMedId, itemId);
//            try {
//                store(pub);
//            } catch (ObjectStoreException e) {
//                throw new ObjectStoreException(e);
//            }
//        }
//        return itemId;
//    }



    /**
     * Create and store a organism item on the first time called.
     *
     * @throws ObjectStoreException os
     */
    protected void createOrganismItem() throws ObjectStoreException {
        org = createItem("Organism");
        org.setAttribute("taxonId", ATH_TAXID);
        store(org);
    }



}
