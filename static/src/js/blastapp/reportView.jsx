import React, { useState, useEffect } from 'react';
import { makeStyles } from '@material-ui/core/styles';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import DropdownButton from 'react-bootstrap/DropdownButton';
import Dropdown from 'react-bootstrap/Dropdown';
import Alert from 'react-bootstrap/Alert';
import BlastReportDescription from './reportDescription.jsx';
import BlastReportAlignments from './reportAlignments.jsx';

const useStyles = makeStyles(theme => ({
    root: {
        width: '100%',
        marginTop: '20px',
        overflowX: 'auto',
        position: 'sticky',
        top: '0px',
        'zIndex': '1000',
    },
}));

function BlastReport({ report }) {
    const [tab, setTab] = useState(0);
    const [algnItem, setAlgnItem] = useState(null);

    function handleTabClick(event, newValue) {
        setTab(newValue);
    }

    useEffect(() => {
        if (algnItem !== null && tab === 1) {
            $('html, body').animate(
                {
                    scrollTop: $('#dln_' + algnItem).offset().top - 60,
                },
                300
            );
            setAlgnItem(null);
        }
    }, [tab]);

    function handleItemClick(event) {
        event.preventDefault();
        setTab(1);
        setAlgnItem(event.target.getAttribute('href'));
    }

    const classes = useStyles();

    if (report.report.results.search.hits.length < 1) {
        return (
            <div className="mt-4 text-align-center">
                <Alert dismissible variant="info">
                    <Alert.Heading>
                        There were no hits for this query...
                    </Alert.Heading>
                    <p>
                        You may try again with a new sequence. Please confirm
                        that you are entering either amino acid sequence(s) or
                        nucleotide sequence(s), (but not both in the same FASTA
                        entry).
                    </p>
                </Alert>
            </div>
        );
    }

    return (
        <div>
            <Paper square className={classes.root}>
                <Tabs
                    value={tab}
                    onChange={handleTabClick}
                    indicatorColor="primary"
                    textColor="primary"
                >
                    <Tab label="Descriptions" />
                    <Tab label="Alignments" />
                </Tabs>
            </Paper>
            {tab === 0 && (
                <div>
                    <BlastReportDescription
                        report={report.report.results}
                        onClick={handleItemClick}
                    />
                </div>
            )}
            {tab === 1 && (
                <div>
                    <BlastReportAlignments report={report.report.results} />
                </div>
            )}
        </div>
    );
}

export default BlastReport;
