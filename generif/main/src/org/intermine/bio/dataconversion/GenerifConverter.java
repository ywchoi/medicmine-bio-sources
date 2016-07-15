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
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.apache.tools.ant.BuildException;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.metadata.StringUtil;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.FormattedTextParser;
import org.intermine.xml.full.Item;

/**
 *
 * @author sc
 */
public class GenerifConverter extends BioFileConverter
{
    private static final Logger LOG = Logger.getLogger(GenerifConverter.class);

    private static final String DATASET_TITLE = "GeneRIF";
    private static final String DATA_SOURCE_NAME = "NCBI";

    private IdResolver rslv;
    private Set<String> taxonIds = new HashSet<String>();

    private Map<String, String> orgItems = new HashMap<String, String>();
    private Map<String, String> pubItems = new HashMap<String, String>();
    private Map<String, String> geneItems = new HashMap<String, String>();

    /**
     * Constructor
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     * @throws ObjectStoreException if problems..
     */
    public GenerifConverter(ItemWriter writer, Model model)
        throws ObjectStoreException {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * Sets the list of taxonIds that should be processed.  All genes will be loaded.
     *
     * @param taxonIds a space-separated list of taxonIds
     */
    public void setGenerifOrganisms(String taxonIds) throws ObjectStoreException {
        this.taxonIds = new HashSet<String>(Arrays.asList(StringUtil.split(taxonIds, " ")));
        LOG.info("Setting list of organisms to " + taxonIds);
        for (String taxonId : this.taxonIds) {
            createOrganismItem(taxonId);
        }
    }

    /**
     *
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
        File currentFile = getCurrentFile();

        if ("generifs_basic".equals(currentFile.getName())) {
            processFile(reader);
        } else {
            LOG.info("WWSS skipping file: " + currentFile.getName());
            //            throw new IllegalArgumentException("Unexpected file: "
            //          + currentFile.getName());
        }
    }

    /**
     * Process all rows of the generifs_basic file, available at
     * ftp://ftp.ncbi.nih.gov/gene/GeneRIF/generifs_basic
     *
     * @param reader
     *            a reader for the generifs_basic file
     *
     * @throws IOException
     * @throws ObjectStoreException
     *
     * FILE FORMAT tsv
     * HEADER
     * #TaxID GeneID PubMedID lastUpdateTimeStamp   GeneRIFText
     * EXAMPLE
     * 3702    814572  17550895        2010-01-21 00:00        AtCCMA interacts with AtCcmB [..]
     *
     *
     */
    private void processFile(Reader reader)
        throws IOException, ObjectStoreException {
        Iterator<?> tsvIter;
        try {
            tsvIter = FormattedTextParser.parseTabDelimitedReader(reader);
        } catch (Exception e) {
            throw new BuildException("cannot parse file: " + getCurrentFile(), e);
        }

        if (taxonIds.isEmpty()) {
            LOG.warn("generif.organisms property not set in project XML file");
        }

        //Create id resolver
        if (rslv == null) {
            rslv = IdResolverService.getIdResolverByOrganism(taxonIds);
        }

        String pid = null;

        int lineNumber = 0;

        while (tsvIter.hasNext()) {
            String[] line = (String[]) tsvIter.next();

            // this can be omitted
            if (lineNumber == 0) {
                checkHeader(line);
                lineNumber++;
                continue;
            }

            String taxid = line[0];
            String geneId = line [1];
            String pubMedId = line [2];
            String timeStamp = line[3];
            String annotation = line [4];

            if (!taxonIds.contains(taxid)) {
                continue;
            }

            int resCount = rslv.countResolutions(taxid, geneId);
            if (resCount != 1) {
                LOG.info("RESOLVER: failed to resolve gene to one identifier, ignoring gene: "
                        + geneId + " count: " + resCount);
                continue;
            }

            // NOT WORKING!? see
            // http://intermine.readthedocs.org/en/latest/database/data-sources/id-resolvers
            //          if (pid == null) {
            //              LOG.info("MISSING ID: " + geneId);
            //              continue;
            //          }


            pid = rslv.resolveId(taxid, geneId).iterator().next();
            LOG.info("READING " + taxid + ": " + pid + "<->" + geneId + "|" + pubMedId + "|"
                    + timeStamp + "--" + annotation);

            Item ann = createGeneRIF(annotation, timeStamp);
            createBioEntity(pid, "Gene");
            createPublication(pubMedId);
            ann.setReference("gene", geneItems.get(pid));
            ann.setReference("organism", orgItems.get(taxid));
            ann.setReference("publication", pubItems.get(pubMedId));
            store(ann);

            lineNumber++;
        }
    }

    /**
     * @param line
     */
    private void checkHeader(String[] line) {
        // column headers - strip off any extra columns - FlyAtlas
        // not necessary
        String[] headers;
        int end = 0;
        for (int i = 0; i < line.length; i++) {
            // if (StringUtils.isEmpty(line[i])) {
            if (line[i].isEmpty()) {
                break;
            }
            end++;
        }
        headers = new String[end];
        System.arraycopy(line, 0, headers, 0, end);
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
     * @param type the type of bioentity (gene, exon..)
     * @throws ObjectStoreException
     */
    private void createBioEntity(String primaryId, String type) throws ObjectStoreException {
        // doing only genes here
        Item bioentity = null;

        if ("Gene".equals(type)) {
            if (!geneItems.containsKey(primaryId)) {
                bioentity = createItem("Gene");
                bioentity.setAttribute("primaryIdentifier", primaryId);
                store(bioentity);
                geneItems.put(primaryId, bioentity.getIdentifier());
            }
        }
    }
    /**
     * Create and store a Publication item on the first time called.
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


    /**
     * Create and store a organism item on the first time called.
     *
     * @param taxonId the taxon ID
     * @throws ObjectStoreException os
     */
    private void createOrganismItem(String taxonId) throws ObjectStoreException {
        if (!orgItems.containsKey(taxonId)) {
            Item org = createItem("Organism");
            org.setAttribute("taxonId", taxonId);
            store(org);
            orgItems.put(taxonId, org.getIdentifier());
        }
    }

}
